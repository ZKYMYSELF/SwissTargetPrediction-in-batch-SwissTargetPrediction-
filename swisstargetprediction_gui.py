#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import re
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from queue import Empty, Queue
from tkinter import END, DISABLED, NORMAL, DoubleVar, StringVar, Tk, Text, filedialog, messagebox, ttk

import swisstargetprediction_batch as core


def split_cells(text: str) -> list[str]:
    items: list[str] = []
    for raw_line in text.replace("\r", "\n").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = [part.strip() for part in re.split(r"\t+", line) if part.strip()]
        if parts:
            items.append(parts[0])
    return items


def looks_like_name_header(value: str) -> bool:
    return value.strip().lower() in {
        "name",
        "compound",
        "compound name",
        "compound_name",
        "drug",
        "化合物名称",
        "名称",
        "物质名称",
    }


def looks_like_smiles_header(value: str) -> bool:
    return value.strip().lower() in {
        "smiles",
        "smi",
        "smile",
        "smiles号",
    }


def parse_column_text(raw: str, kind: str) -> list[str]:
    items = split_cells(raw)
    if len(items) > 1:
        if kind == "name" and looks_like_name_header(items[0]):
            items = items[1:]
        if kind == "smiles" and looks_like_smiles_header(items[0]):
            items = items[1:]
    return items


def parse_clipboard_two_columns(text: str) -> tuple[list[str], list[str]]:
    rows: list[list[str]] = []
    for raw_line in text.replace("\r", "\n").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = [part.strip() for part in re.split(r"\t+", line) if part.strip()]
        if parts:
            rows.append(parts)

    if not rows:
        return [], []

    if len(rows[0]) >= 2 and looks_like_name_header(rows[0][0]) and looks_like_smiles_header(rows[0][1]):
        rows = rows[1:]

    names: list[str] = []
    smiles: list[str] = []
    for row in rows:
        if len(row) < 2:
            continue
        names.append(row[0])
        smiles.append(row[1])
    return names, smiles


def sanitize_filename(value: str, max_length: int = 60) -> str:
    value = re.sub(r'[<>:"/\\|?*]+', "_", value).strip()
    value = re.sub(r"\s+", " ", value).strip(" .")
    if not value:
        value = "未命名"
    if len(value) > max_length:
        value = value[:max_length].rstrip(" ._")
    return value


def timestamp_text() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def build_output_rows(
    compound_name: str,
    sequence_number: int,
    smiles: str,
    organism: str,
    result_url: str,
    prediction_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in prediction_rows:
        rows.append(
            {
                "化合物名称": compound_name,
                "序号": str(sequence_number),
                "SMILES": smiles,
                "物种": organism,
                "target": row.get("target", ""),
                "common_name": row.get("common_name", ""),
                "uniprot_id": row.get("uniprot_id", ""),
                "chembl_id": row.get("chembl_id", ""),
                "target_class": row.get("target_class", ""),
                "probability": row.get("probability", ""),
                "known_actives_3d": row.get("known_actives_3d", ""),
                "known_actives_2d": row.get("known_actives_2d", ""),
                "target_link": row.get("target_link", ""),
                "common_name_link": row.get("common_name_link", ""),
                "uniprot_link": row.get("uniprot_link", ""),
                "chembl_link": row.get("chembl_link", ""),
                "known_actives_link": row.get("known_actives_link", ""),
                "result_url": result_url,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise RuntimeError("没有可写入的预测结果。")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def extract_site_error_message(page_text: str) -> str:
    match = re.search(r"(?is)<div[^>]*id=['\"]content['\"][^>]*>(.*?)</div>", page_text)
    snippet = match.group(1) if match else page_text
    cleaned = re.sub(r"(?is)<script.*?</script>", " ", snippet)
    cleaned = re.sub(r"(?is)<style.*?</style>", " ", cleaned)
    cleaned = re.sub(r"(?is)<[^>]+>", " ", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return cleaned[:500] if cleaned else "SwissTargetPrediction 返回了错误页，但没有可读提示。"


@dataclass
class JobPair:
    index: int
    name: str
    smiles: str


@dataclass
class JobResult:
    index: int
    name: str
    smiles: str
    status: str
    message: str
    file_name: str = ""


class SwissTargetPredictionGUI:
    def __init__(self, root: Tk) -> None:
        self.root = root
        self.root.title("SwissTargetPrediction 批量靶点预测")
        self.root.geometry("1600x1000")
        self.root.minsize(1280, 860)

        self.message_queue: Queue[tuple[str, object]] = Queue()
        self.worker_thread: threading.Thread | None = None
        self.running = False
        self.hover_scroll_widget: object | None = None
        self.app_dir = Path(__file__).resolve().parent

        self.organism_var = StringVar(value="Homo_sapiens")
        self.output_dir_var = StringVar(value=str(self.app_dir / "输出结果"))
        self.delay_var = DoubleVar(value=8.0)
        self.status_var = StringVar(value="请输入化合物名称和 SMILES。")
        self.count_var = StringVar(value="名称 0 条，SMILES 0 条")
        self.validation_var = StringVar(value="等待输入。")
        self.summary_var = StringVar(value="完成 0 个，失败 0 个")

        self._preview_after_id: str | None = None
        self.row_lookup: dict[int, str] = {}

        self._build_style()
        self._build_ui()
        self._poll_queue()

    def _build_style(self) -> None:
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure("Title.TLabel", font=("Microsoft YaHei UI", 16, "bold"))
        style.configure("Section.TLabelframe.Label", font=("Microsoft YaHei UI", 10, "bold"))

    def _build_ui(self) -> None:
        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True)

        from tkinter import Canvas

        self.page_canvas = Canvas(outer, highlightthickness=0)
        self.page_scrollbar = ttk.Scrollbar(outer, orient="vertical", command=self.page_canvas.yview)
        self.page_canvas.configure(yscrollcommand=self.page_scrollbar.set)

        self.page_scrollbar.pack(side="right", fill="y")
        self.page_canvas.pack(side="left", fill="both", expand=True)

        self.page_content = ttk.Frame(self.page_canvas, padding=14)
        self.page_window = self.page_canvas.create_window((0, 0), window=self.page_content, anchor="nw")

        self.page_content.bind("<Configure>", self._on_page_content_configure)
        self.page_canvas.bind("<Configure>", self._on_page_canvas_configure)
        self.root.bind_all("<MouseWheel>", self._on_mousewheel, add="+")
        self.root.bind_all("<Button-4>", self._on_linux_scroll, add="+")
        self.root.bind_all("<Button-5>", self._on_linux_scroll, add="+")

        ttk.Label(self.page_content, text="SwissTargetPrediction 批量靶点预测", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            self.page_content,
            text="左边输入化合物名称，右边输入 SMILES。支持直接从 Excel 复制两列后粘贴导入。",
            foreground="#555555",
        ).pack(anchor="w", pady=(4, 10))

        control_frame = ttk.LabelFrame(self.page_content, text="全局设置", style="Section.TLabelframe", padding=10)
        control_frame.pack(fill="x", pady=(0, 10))

        ttk.Label(control_frame, text="物种").grid(row=0, column=0, sticky="w", padx=(0, 6))
        ttk.Combobox(
            control_frame,
            textvariable=self.organism_var,
            values=["Homo_sapiens", "Mus_musculus", "Rattus_norvegicus"],
            state="readonly",
            width=24,
        ).grid(row=0, column=1, sticky="w", padx=(0, 12))

        ttk.Label(control_frame, text="提交间隔(秒)").grid(row=0, column=2, sticky="w", padx=(0, 6))
        ttk.Spinbox(control_frame, from_=0.0, to=60.0, increment=0.5, textvariable=self.delay_var, width=10).grid(
            row=0, column=3, sticky="w", padx=(0, 12)
        )

        ttk.Label(control_frame, text="输出目录").grid(row=0, column=4, sticky="w", padx=(0, 6))
        ttk.Entry(control_frame, textvariable=self.output_dir_var, width=64).grid(row=0, column=5, sticky="ew", padx=(0, 8))
        ttk.Button(control_frame, text="选择", command=self.choose_output_dir).grid(row=0, column=6, sticky="w")
        control_frame.columnconfigure(5, weight=1)

        button_row = ttk.Frame(self.page_content)
        button_row.pack(fill="x", pady=(0, 10))
        ttk.Button(button_row, text="从 Excel 剪贴板导入两列", command=self.import_from_clipboard).pack(side="left")
        ttk.Button(button_row, text="清空全部", command=self.clear_all).pack(side="left", padx=8)
        self.start_button = ttk.Button(button_row, text="开始预测并逐个保存", command=self.start_prediction)
        self.start_button.pack(side="left", padx=8)
        ttk.Button(button_row, text="打开输出目录", command=self.open_output_dir).pack(side="left", padx=8)

        input_frame = ttk.Frame(self.page_content)
        input_frame.pack(fill="both", expand=False)

        self._build_column(input_frame, "化合物名称", 0)
        self._build_column(input_frame, "SMILES", 1)

        status_frame = ttk.LabelFrame(self.page_content, text="状态", style="Section.TLabelframe", padding=10)
        status_frame.pack(fill="x", pady=(10, 10))
        ttk.Label(status_frame, textvariable=self.count_var).pack(anchor="w")
        ttk.Label(status_frame, textvariable=self.validation_var, foreground="#aa0000").pack(anchor="w", pady=(4, 0))
        ttk.Label(status_frame, textvariable=self.status_var, foreground="#004080").pack(anchor="w", pady=(4, 0))
        ttk.Label(status_frame, textvariable=self.summary_var, foreground="#006600").pack(anchor="w", pady=(4, 0))

        preview_frame = ttk.LabelFrame(self.page_content, text="配对预览", style="Section.TLabelframe", padding=10)
        preview_frame.pack(fill="both", expand=True, pady=(0, 10))
        self.preview_tree = ttk.Treeview(
            preview_frame,
            columns=("序号", "化合物名称", "SMILES", "状态", "说明", "文件"),
            show="headings",
            height=12,
        )
        for column, width in [
            ("序号", 70),
            ("化合物名称", 240),
            ("SMILES", 420),
            ("状态", 90),
            ("说明", 560),
            ("文件", 260),
        ]:
            self.preview_tree.heading(column, text=column)
            self.preview_tree.column(column, width=width, anchor="w")
        self.preview_tree.tag_configure("ok", foreground="#0a7a0a")
        self.preview_tree.tag_configure("fail", foreground="#b00020")
        self.preview_tree.tag_configure("pending", foreground="#666666")
        preview_scroll_y = ttk.Scrollbar(preview_frame, orient="vertical", command=self.preview_tree.yview)
        preview_scroll_x = ttk.Scrollbar(preview_frame, orient="horizontal", command=self.preview_tree.xview)
        self.preview_tree.configure(yscrollcommand=preview_scroll_y.set, xscrollcommand=preview_scroll_x.set)
        self.preview_tree.grid(row=0, column=0, sticky="nsew")
        preview_scroll_y.grid(row=0, column=1, sticky="ns")
        preview_scroll_x.grid(row=1, column=0, sticky="ew")
        preview_frame.rowconfigure(0, weight=1)
        preview_frame.columnconfigure(0, weight=1)
        self._register_scrollable(self.preview_tree)

        log_frame = ttk.LabelFrame(self.page_content, text="运行日志", style="Section.TLabelframe", padding=10)
        log_frame.pack(fill="both", expand=True)
        self.log_text = Text(log_frame, height=10, wrap="word")
        log_scroll = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=log_scroll.set)
        self.log_text.pack(side="left", fill="both", expand=True)
        log_scroll.pack(side="right", fill="y")
        self._register_scrollable(self.log_text)

        self.name_text.focus_set()
        self.refresh_preview()

    def _build_column(self, parent: ttk.Frame, title: str, column_index: int) -> None:
        frame = ttk.LabelFrame(parent, text=title, style="Section.TLabelframe", padding=10)
        frame.grid(row=0, column=column_index, sticky="nsew", padx=(0, 10) if column_index == 0 else (10, 0))
        parent.columnconfigure(column_index, weight=1)

        ttk.Label(frame, text="直接粘贴即可。支持 Excel 复制后的换行/制表符文本。", foreground="#666666").pack(anchor="w")

        text_frame = ttk.Frame(frame)
        text_frame.pack(fill="both", expand=True, pady=(6, 6))
        text_widget = Text(text_frame, height=14, wrap="none", undo=True)
        text_scroll = ttk.Scrollbar(text_frame, orient="vertical", command=text_widget.yview)
        text_widget.configure(yscrollcommand=text_scroll.set)
        text_widget.pack(side="left", fill="both", expand=True)
        text_scroll.pack(side="right", fill="y")
        text_widget.bind("<<Modified>>", self.on_text_changed)
        text_widget.bind("<KeyRelease>", self.on_text_changed)
        text_widget.bind("<FocusOut>", self.on_text_changed)
        self._register_scrollable(text_widget)

        preview = Text(frame, height=8, wrap="none", state=DISABLED, bg="#f7f7f7")
        preview_scroll = ttk.Scrollbar(frame, orient="vertical", command=preview.yview)
        preview.configure(yscrollcommand=preview_scroll.set)
        preview.pack(fill="both", expand=False, side="left")
        preview_scroll.pack(side="right", fill="y")
        self._register_scrollable(preview)

        if column_index == 0:
            self.name_text = text_widget
            self.name_preview = preview
        else:
            self.smiles_text = text_widget
            self.smiles_preview = preview

    def _register_scrollable(self, widget: object) -> None:
        try:
            widget.bind("<Enter>", lambda event, target=widget: self._set_hover_scroll_widget(target), add="+")  # type: ignore[attr-defined]
            widget.bind("<Leave>", lambda event, target=widget: self._clear_hover_scroll_widget(target), add="+")  # type: ignore[attr-defined]
        except Exception:
            pass

    def _set_hover_scroll_widget(self, widget: object) -> None:
        self.hover_scroll_widget = widget

    def _clear_hover_scroll_widget(self, widget: object) -> None:
        if self.hover_scroll_widget is widget:
            self.hover_scroll_widget = None

    def _event_delta(self, event) -> int:
        if getattr(event, "num", None) == 4:
            return 1
        if getattr(event, "num", None) == 5:
            return -1
        if getattr(event, "delta", 0) == 0:
            return 0
        return int(event.delta / 120)

    def _can_scroll(self, widget: object, delta: int) -> bool:
        try:
            first, last = widget.yview()  # type: ignore[attr-defined]
            if delta > 0 and first <= 0.0:
                return False
            if delta < 0 and last >= 1.0:
                return False
            return True
        except Exception:
            return False

    def _scroll_widget(self, widget: object, delta: int) -> bool:
        if not self._can_scroll(widget, delta):
            return False
        try:
            widget.yview_scroll(-delta, "units")  # type: ignore[attr-defined]
            return True
        except Exception:
            return False

    def _on_mousewheel(self, event) -> str | None:
        delta = self._event_delta(event)
        if delta == 0:
            return None
        if self.hover_scroll_widget is not None and self._scroll_widget(self.hover_scroll_widget, delta):
            return "break"
        if self._scroll_widget(self.page_canvas, delta):
            return "break"
        return None

    def _on_linux_scroll(self, event) -> str | None:
        return self._on_mousewheel(event)

    def _on_page_content_configure(self, event) -> None:
        self.page_canvas.configure(scrollregion=self.page_canvas.bbox("all"))

    def _on_page_canvas_configure(self, event) -> None:
        self.page_canvas.itemconfigure(self.page_window, width=event.width)

    def choose_output_dir(self) -> None:
        selected = filedialog.askdirectory(initialdir=self.output_dir_var.get() or str(self.app_dir))
        if selected:
            self.output_dir_var.set(selected)

    def open_output_dir(self) -> None:
        path = Path(self.output_dir_var.get())
        path.mkdir(parents=True, exist_ok=True)
        try:
            os.startfile(path)
        except Exception as exc:
            messagebox.showerror("打开失败", str(exc))

    def clear_all(self) -> None:
        self.name_text.delete("1.0", END)
        self.smiles_text.delete("1.0", END)
        self.refresh_preview()

    def import_from_clipboard(self) -> None:
        try:
            clipboard_text = self.root.clipboard_get()
        except Exception:
            messagebox.showerror(
                "剪贴板导入失败",
                "剪贴板里没有可读取的文本。\n"
                "请先在 Excel 里复制两列文本，再点这个按钮。\n"
                "如果你复制的是图片/截图，这个版本暂时不能直接识别。"
            )
            return

        names, smiles = parse_clipboard_two_columns(clipboard_text)
        if not names or not smiles or len(names) != len(smiles):
            messagebox.showerror(
                "剪贴板格式不对",
                "请从 Excel 复制两列数据：左列化合物名称，右列 SMILES。\n"
                "每一行必须是一对。"
            )
            return

        self.name_text.delete("1.0", END)
        self.name_text.insert("1.0", "\n".join(names))
        self.smiles_text.delete("1.0", END)
        self.smiles_text.insert("1.0", "\n".join(smiles))
        self.refresh_preview()
        self.log(f"已从剪贴板导入 {len(names)} 对数据。")

    def on_text_changed(self, event=None) -> None:
        widget = getattr(event, "widget", None)
        if widget is not None and isinstance(widget, Text):
            widget.edit_modified(False)
        if self._preview_after_id is not None:
            self.root.after_cancel(self._preview_after_id)
        self._preview_after_id = self.root.after(120, self.refresh_preview)

    def get_name_items(self) -> list[str]:
        return parse_column_text(self.name_text.get("1.0", END), "name")

    def get_smiles_items(self) -> list[str]:
        return parse_column_text(self.smiles_text.get("1.0", END), "smiles")

    def refresh_preview(self) -> None:
        self._preview_after_id = None
        names = self.get_name_items()
        smiles = self.get_smiles_items()
        self.count_var.set(f"名称 {len(names)} 条，SMILES {len(smiles)} 条")

        if not names and not smiles:
            self.validation_var.set("等待输入。")
        elif len(names) != len(smiles):
            self.validation_var.set("数量不一致，无法匹配。请检查两边条目数。")
        else:
            self.validation_var.set("数量匹配，可以开始预测。")

        self._render_side_preview(self.name_preview, names, "化合物名称")
        self._render_side_preview(self.smiles_preview, smiles, "SMILES")

        self.preview_tree.delete(*self.preview_tree.get_children())
        self.row_lookup.clear()
        for index, name in enumerate(names, start=1):
            smile = smiles[index - 1] if index - 1 < len(smiles) else ""
            iid = str(index)
            self.row_lookup[index] = iid
            self.preview_tree.insert(
                "",
                END,
                iid=iid,
                values=(index, name, smile, "待处理", "", ""),
                tags=("pending",),
            )

        enabled = bool(names and smiles and len(names) == len(smiles) and not self.running)
        self.start_button.configure(state=NORMAL if enabled else DISABLED)

    def _render_side_preview(self, widget: Text, items: list[str], label: str) -> None:
        widget.configure(state=NORMAL)
        widget.delete("1.0", END)
        if items:
            widget.insert("1.0", "\n".join(f"{index}. {item}" for index, item in enumerate(items, start=1)))
        else:
            widget.insert("1.0", f"这里会显示编号后的{label}预览。")
        widget.configure(state=DISABLED)

    def log(self, message: str) -> None:
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.configure(state=NORMAL)
        self.log_text.insert(END, f"[{timestamp}] {message}\n")
        self.log_text.see(END)
        self.log_text.configure(state=DISABLED)

    def update_preview_row(self, result: JobResult) -> None:
        iid = self.row_lookup.get(result.index)
        if not iid:
            return
        self.preview_tree.item(
            iid,
            values=(result.index, result.name, result.smiles, result.status, result.message, result.file_name),
            tags=("ok" if result.status == "成功" else "fail" if result.status == "失败" else "pending",),
        )

    def start_prediction(self) -> None:
        if self.running:
            return

        names = self.get_name_items()
        smiles = self.get_smiles_items()
        if len(names) != len(smiles):
            messagebox.showerror("数量不一致", "名称和 SMILES 的数量必须一致。")
            return
        if not names:
            messagebox.showerror("没有数据", "请先输入至少一对化合物名称和 SMILES。")
            return

        out_dir = Path(self.output_dir_var.get())
        out_dir.mkdir(parents=True, exist_ok=True)

        self.running = True
        self.start_button.configure(state=DISABLED)
        self.status_var.set("开始预测，正在逐个提交。")
        self.summary_var.set(f"完成 0 个，失败 0 个")
        self.log(f"开始任务，共 {len(names)} 对，物种 {self.organism_var.get()}，输出到 {out_dir}.")

        pairs = [JobPair(index=i, name=name, smiles=smile) for i, (name, smile) in enumerate(zip(names, smiles), start=1)]
        self.worker_thread = threading.Thread(target=self.worker_run, args=(pairs, out_dir), daemon=True)
        self.worker_thread.start()

    def submit_query_with_retry(
        self,
        session,
        pair: JobPair,
        organism: str,
        max_submit_attempts: int = 3,
    ) -> str:
        submit_error: Exception | None = None
        for attempt in range(1, max_submit_attempts + 1):
            try:
                return core.submit_query(session, core.Query(smiles=pair.smiles, name=pair.name), organism)
            except Exception as exc:
                submit_error = exc
                if attempt < max_submit_attempts:
                    backoff_seconds = 5 * attempt
                    self.message_queue.put((
                        "log", f"{pair.index} {pair.name}: ?????? {attempt} ????{backoff_seconds} ?????"
                    ))
                    time.sleep(backoff_seconds)
        raise RuntimeError(
            f"{submit_error}???? {max_submit_attempts} ??????????????? SwissTargetPrediction ???????"
        )

    def worker_run(
        self, pairs: list[JobPair], out_dir: Path) -> None:
        session = core.make_session()
        organism = self.organism_var.get()
        delay_seconds = float(self.delay_var.get())
        run_timestamp = timestamp_text()
        success_count = 0
        failure_count = 0
        max_submit_attempts = 3

        for pair in pairs:
            self.message_queue.put(("status", f"正在预测 {pair.index}/{len(pairs)}: {pair.name}"))
            try:
                result_url = self.submit_query_with_retry(session, pair, organism, max_submit_attempts)
                self.message_queue.put(("log", f"{pair.index} {pair.name}: 已提交，等待结果。"))
                page_html = core.wait_for_result(session, result_url, timeout_seconds=180.0, wait_seconds=10.0)
                prediction_rows = core.parse_prediction_rows(page_html)
                output_rows = build_output_rows(pair.name, pair.index, pair.smiles, organism, result_url, prediction_rows)
                file_name = f"{pair.index}_{sanitize_filename(pair.name)}_{run_timestamp}.csv"
                output_path = out_dir / file_name
                write_csv(output_path, output_rows)
                success_count += 1
                self.message_queue.put(("row", JobResult(pair.index, pair.name, pair.smiles, "成功", f"已保存 {len(output_rows)} 行", output_path.name)))
                self.message_queue.put(("log", f"{pair.index} {pair.name}: 已保存 {len(output_rows)} 行 -> {output_path.name}"))
            except Exception as exc:
                failure_count += 1
                message = str(exc)
                error_url_match = re.search(r"https?://\S*error_page\.php\?\S+", message)
                if error_url_match:
                    try:
                        error_url = error_url_match.group(0).rstrip(").,;")
                        page_text = session.get(error_url, timeout=60).text
                        site_message = extract_site_error_message(page_text)
                        message = f"{message} | 网站提示: {site_message}"
                    except Exception:
                        pass
                self.message_queue.put(("row", JobResult(pair.index, pair.name, pair.smiles, "失败", message, "")))
                self.message_queue.put(("log", f"{pair.index} {pair.name}: 失败 - {message}"))

            if pair.index < len(pairs) and delay_seconds > 0:
                time.sleep(delay_seconds)

        self.message_queue.put(("finished", (success_count, failure_count, str(out_dir))))

    def _poll_queue(self) -> None:
        try:
            while True:
                kind, payload = self.message_queue.get_nowait()
                if kind == "log":
                    self.log(str(payload))
                elif kind == "status":
                    self.status_var.set(str(payload))
                elif kind == "row":
                    self.update_preview_row(payload)  # type: ignore[arg-type]
                elif kind == "finished":
                    success_count, failure_count, out_dir = payload  # type: ignore[misc]
                    self.running = False
                    self.status_var.set(f"完成：成功 {success_count} 个，失败 {failure_count} 个。")
                    self.summary_var.set(f"完成 {success_count} 个，失败 {failure_count} 个")
                    self.log(f"全部完成：成功 {success_count} 个，失败 {failure_count} 个。")
                    messagebox.showinfo("完成", f"预测结束。\n成功 {success_count} 个，失败 {failure_count} 个。\n输出目录：{out_dir}")
                    enabled = bool(self.get_name_items() and self.get_smiles_items() and len(self.get_name_items()) == len(self.get_smiles_items()))
                    self.start_button.configure(state=NORMAL if enabled else DISABLED)
                self.message_queue.task_done()
        except Empty:
            pass
        self.root.after(120, self._poll_queue)


def main() -> int:
    root = Tk()
    SwissTargetPredictionGUI(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
