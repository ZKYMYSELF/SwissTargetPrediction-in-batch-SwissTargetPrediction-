#!/usr/bin/env python3
"""Batch SwissTargetPrediction client.

Usage examples:
  python swisstargetprediction_batch.py --smiles "CCC[C@H]1CCCCN1"
  python swisstargetprediction_batch.py --input compounds.txt --output results.csv
  python swisstargetprediction_batch.py --input compounds.csv --smiles-column smiles --name-column name

The script submits SMILES to SwissTargetPrediction, follows the job redirect,
waits for the results table to appear, and exports all predicted targets to CSV.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import html
import os
import re
import sys
import time
from pathlib import Path
from typing import Iterable, Iterator

import requests


BASE_URL = "https://www.swisstargetprediction.ch"
SUBMIT_URL = f"{BASE_URL}/predict.php"
INDEX_URL = f"{BASE_URL}/index.php"

DEFAULT_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
}


@dataclasses.dataclass
class Query:
    smiles: str
    name: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch SwissTargetPrediction scraper for SMILES strings."
    )
    parser.add_argument("--smiles", help="Single SMILES to submit.")
    parser.add_argument(
        "--name",
        default="",
        help="Optional name for the single SMILES query.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Input file. Plain text: one SMILES per line. CSV: use --smiles-column.",
    )
    parser.add_argument(
        "--smiles-column",
        default="smiles",
        help="CSV column containing SMILES (case-insensitive).",
    )
    parser.add_argument(
        "--name-column",
        default="",
        help="Optional CSV column containing a query name (case-insensitive).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("swisstargetprediction_results.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--organism",
        choices=["Homo_sapiens", "Mus_musculus", "Rattus_norvegicus"],
        default="Homo_sapiens",
        help="Target species.",
    )
    parser.add_argument(
        "--wait-seconds",
        type=float,
        default=10.0,
        help="Polling interval while waiting for the result page.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=180.0,
        help="Maximum time to wait for one prediction job.",
    )
    return parser.parse_args()


def parse_args_or_interactive() -> argparse.Namespace:
    if len(sys.argv) > 1:
        return parse_args()

    print("SwissTargetPrediction 批量工具")
    print("1) 单个 SMILES")
    print("2) 批量文件(txt/csv)")
    choice = input("请选择 1 或 2: ").strip()

    if choice == "1":
        smiles = input("请输入 SMILES: ").strip()
        if not smiles:
            raise SystemExit("SMILES 不能为空。")
        name = input("可选名称(直接回车跳过): ").strip()
        output = input("输出 CSV 文件名(默认 swisstargetprediction_results.csv): ").strip()
        namespace = argparse.Namespace(
            smiles=smiles,
            name=name,
            input=None,
            smiles_column="smiles",
            name_column="",
            output=Path(output) if output else Path("swisstargetprediction_results.csv"),
            organism="Homo_sapiens",
            wait_seconds=10.0,
            timeout_seconds=180.0,
        )
        return namespace

    if choice == "2":
        input_path = input("请输入 txt/csv 文件路径: ").strip().strip('"')
        if not input_path:
            raise SystemExit("文件路径不能为空。")
        output = input("输出 CSV 文件名(默认 swisstargetprediction_results.csv): ").strip()
        organism = input("物种 [1=Homo_sapiens, 2=Mus_musculus, 3=Rattus_norvegicus] (默认 1): ").strip()
        organism_map = {
            "1": "Homo_sapiens",
            "2": "Mus_musculus",
            "3": "Rattus_norvegicus",
        }
        namespace = argparse.Namespace(
            smiles=None,
            name="",
            input=Path(input_path),
            smiles_column="smiles",
            name_column="",
            output=Path(output) if output else Path("swisstargetprediction_results.csv"),
            organism=organism_map.get(organism, "Homo_sapiens"),
            wait_seconds=10.0,
            timeout_seconds=180.0,
        )
        return namespace

    raise SystemExit("已取消。")


def load_queries(args: argparse.Namespace) -> list[Query]:
    if args.smiles:
        return [Query(smiles=args.smiles.strip(), name=args.name.strip())]
    if not args.input:
        raise SystemExit("Provide either --smiles or --input.")
    if not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")

    suffix = args.input.suffix.lower()
    if suffix == ".csv":
        return load_queries_from_csv(args.input, args.smiles_column, args.name_column)
    return load_queries_from_text(args.input)


def load_queries_from_text(path: Path) -> list[Query]:
    queries: list[Query] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        queries.append(Query(smiles=line))
    return queries


def load_queries_from_csv(path: Path, smiles_column: str, name_column: str) -> list[Query]:
    text = path.read_text(encoding="utf-8-sig")
    sample = text[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error:
        dialect = csv.get_dialect("excel")

    rows = list(csv.reader(text.splitlines(), dialect))
    if not rows:
        return []

    header = [cell.strip().lower() for cell in rows[0]]
    smiles_index = None
    name_index = None
    if smiles_column.lower() in header:
        smiles_index = header.index(smiles_column.lower())
    elif len(rows[0]) >= 1:
        smiles_index = 0

    if name_column and name_column.lower() in header:
        name_index = header.index(name_column.lower())

    queries: list[Query] = []
    data_rows = rows[1:] if any(col in header for col in (smiles_column.lower(), name_column.lower())) else rows
    for row in data_rows:
        if not row or all(not cell.strip() for cell in row):
            continue
        if smiles_index is None or smiles_index >= len(row):
            continue
        smiles = row[smiles_index].strip()
        if not smiles:
            continue
        name = ""
        if name_index is not None and name_index < len(row):
            name = row[name_index].strip()
        queries.append(Query(smiles=smiles, name=name))
    return queries


def make_session() -> requests.Session:
    session = requests.Session()
    session.headers.update(DEFAULT_HEADERS)
    return session


def submit_query(session: requests.Session, query: Query, organism: str) -> str:
    session.get(INDEX_URL, timeout=60)
    response = session.post(
        SUBMIT_URL,
        headers={
            "Referer": INDEX_URL,
            "Content-Type": "application/x-www-form-urlencoded",
        },
        data={
            "smiles": query.smiles,
            "organism": organism,
            "ioi": "2",
        },
        timeout=180,
    )
    redirect = extract_js_redirect(response.text)
    if not redirect:
        raise RuntimeError("SwissTargetPrediction did not return a job redirect.")
    resolved = requests.compat.urljoin(BASE_URL + "/", redirect)
    if "error_page.php" in resolved:
        raise RuntimeError(f"SwissTargetPrediction rejected the submission: {resolved}")
    return resolved


def extract_js_redirect(html_text: str) -> str | None:
    match = re.search(r'location\.replace\("([^"]+)"\)', html_text)
    if match:
        return match.group(1)
    match = re.search(r"location\.replace\('([^']+)'\)", html_text)
    if match:
        return match.group(1)
    return None


def wait_for_result(session: requests.Session, result_url: str, timeout_seconds: float, wait_seconds: float) -> str:
    deadline = time.monotonic() + timeout_seconds
    last_text = ""
    while time.monotonic() < deadline:
        response = session.get(result_url, timeout=180)
        last_text = response.text
        if 'id="resultTable"' in last_text or "id='resultTable'" in last_text:
            return last_text
        time.sleep(wait_seconds)
    raise TimeoutError(
        f"Timed out waiting for result table at {result_url}. Last page size: {len(last_text)} bytes."
    )


def parse_prediction_rows(page_html: str) -> list[dict[str, str]]:
    table_html = extract_table_html(page_html, "resultTable")
    if not table_html:
        raise RuntimeError("Result table not found in the returned page.")

    tbody_match = re.search(r"(?is)<tbody[^>]*>(.*?)</tbody>", table_html)
    body_html = tbody_match.group(1) if tbody_match else table_html
    row_htmls = re.findall(r"(?is)<tr[^>]*>(.*?)</tr>", body_html)
    rows: list[dict[str, str]] = []
    for row_html in row_htmls:
        cells = re.findall(r"(?is)<td[^>]*>(.*?)</td>", row_html)
        if len(cells) < 7:
            continue
        target_cell = cells[0]
        common_name_cell = cells[1]
        uniprot_cell = cells[2]
        chembl_cell = cells[3]
        class_cell = cells[4]
        probability_cell = cells[5]
        known_actives_cell = cells[6]
        probability = extract_probability(probability_cell)
        probability_sort = probability_to_float(probability)
        known_3d, known_2d = extract_known_actives(known_actives_cell)
        rows.append(
            {
                "target": clean_cell_text(target_cell),
                "common_name": clean_cell_text(common_name_cell),
                "uniprot_id": clean_cell_text(uniprot_cell),
                "chembl_id": clean_cell_text(chembl_cell),
                "target_class": clean_cell_text(class_cell),
                "probability": probability,
                "known_actives_3d": known_3d,
                "known_actives_2d": known_2d,
                "target_link": first_href(target_cell),
                "common_name_link": first_href(common_name_cell),
                "uniprot_link": first_href(uniprot_cell),
                "chembl_link": first_href(chembl_cell),
                "known_actives_link": first_href(known_actives_cell),
                "_probability_sort": probability_sort,
            }
        )
    rows.sort(key=lambda row: row.get("_probability_sort", 0.0), reverse=True)
    for row in rows:
        row.pop("_probability_sort", None)
    return rows


def extract_table_html(page_html: str, table_id: str) -> str | None:
    pattern = rf'(?is)<table[^>]*id=["\']{re.escape(table_id)}["\'][^>]*>(.*?)</table>'
    match = re.search(pattern, page_html)
    return match.group(0) if match else None


def extract_probability(cell_html: str) -> str:
    match = re.search(r"(?is)<span[^>]*>\s*([0-9]+(?:\.[0-9]+)?)\s*</span>", cell_html)
    if match:
        return match.group(1)
    matches = re.findall(r"([0-9]+(?:\.[0-9]+)?)", cell_html)
    if matches:
        matches.sort(key=len, reverse=True)
        return matches[0]
    return clean_cell_text(cell_html)


def probability_to_float(probability: str) -> float:
    try:
        return float(probability)
    except (TypeError, ValueError):
        return 0.0


def extract_known_actives(cell_html: str) -> tuple[str, str]:
    match_3d = re.search(r'method=3D[^>]*>(\d+)', cell_html, re.I)
    match_2d = re.search(r'method=FP2[^>]*>(\d+)', cell_html, re.I)
    if match_3d or match_2d:
        return (match_3d.group(1) if match_3d else "", match_2d.group(1) if match_2d else "")

    text = html.unescape(re.sub(r"(?is)<[^>]+>", " ", cell_html))
    values = re.findall(r"\b\d+\b", text)
    if len(values) >= 2:
        return values[0], values[1]
    if len(values) == 1:
        return values[0], ""
    return "", ""


def clean_cell_text(cell_html: str) -> str:
    text = html.unescape(re.sub(r"(?is)<[^>]+>", " ", cell_html))
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def first_href(cell_html: str) -> str:
    match = re.search(r'(?is)href=["\']([^"\']+)["\']', cell_html)
    if not match:
        return ""
    return html.unescape(match.group(1))


def write_output(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise RuntimeError("No prediction rows were parsed.")
    fieldnames = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args_or_interactive()
    queries = load_queries(args)
    if not queries:
        raise SystemExit("No SMILES found in the input.")

    session = make_session()
    all_rows: list[dict[str, str]] = []
    for index, query in enumerate(queries, start=1):
        print(f"[{index}/{len(queries)}] Submitting {query.name or query.smiles}", file=sys.stderr)
        result_url = submit_query(session, query, args.organism)
        print(f"  result: {result_url}", file=sys.stderr)
        page_html = wait_for_result(session, result_url, args.timeout_seconds, args.wait_seconds)
        rows = parse_prediction_rows(page_html)
        for row in rows:
            row["query_name"] = query.name
            row["query_smiles"] = query.smiles
            row["organism"] = args.organism
            row["result_url"] = result_url
        all_rows.extend(rows)
        print(f"  parsed {len(rows)} targets", file=sys.stderr)

    write_output(args.output, all_rows)
    print(f"Wrote {len(all_rows)} rows to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        if len(sys.argv) == 1 and sys.stdin.isatty():
            try:
                input("按回车退出...")
            except EOFError:
                pass
        raise SystemExit(1)
    raise SystemExit(exit_code)
