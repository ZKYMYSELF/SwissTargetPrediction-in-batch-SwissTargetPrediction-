# SwissTargetPrediction Windows Release

This repository provides a Windows-friendly GUI for batch SwissTargetPrediction jobs.

## What it does
- Paste compound names on the left and `SMILES` on the right
- Validate that both sides have the same number of rows
- Import tab-separated text copied from Excel
- Submit compounds one by one to reduce the chance of rate limiting
- Save one CSV per compound with the compound name in column A

## Files
- `gui.py`: thin entry point for the GUI
- `swisstargetprediction_gui.py`: GUI application
- `swisstargetprediction_batch.py`: SwissTargetPrediction submission and parsing logic
- `run_gui.bat`: double-click launcher for Windows
- `requirements.txt`: Python dependency list

## Requirements
- Python 3.9+
- `requests`

Install dependencies:

```bash
pip install -r requirements.txt
```

## Run
- Double-click `run_gui.bat`
- Or run:

```bash
python gui.py
```

## Output
- Each compound is exported to its own CSV file
- Filenames use the format `index_compoundname_timestamp.csv`
- Output is written to the app folder by default, inside `输出结果`

## Notes
- If a compound fails, the app skips that file and continues with the next one
- The preview table marks success and failure separately
- The app tries to show a readable site error message when SwissTargetPrediction returns one

## Credits
- Prepared with Codex
