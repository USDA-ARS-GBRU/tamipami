#!/usr/bin/env python3
"""Summarize everything in results/ into tables and figures, and refresh the README.

Versions are discovered from the ``results/cli_<label>.csv`` / ``app_<label>.csv`` files and their
``.meta.json`` sidecars, so nothing here names a version. Outputs (all in results/):
  report.html (self-contained, presentation tables + figures), cli_summary.md, app_summary.md,
  system.md, cli_benchmark.png, app_benchmark.png
The README block between the AUTO markers is rewritten from those files.
"""
import argparse
import base64
import html
import json
import re
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import bench_common as bc

CLI_METRICS = [("process_time", "Process time (s)"), ("predict_time", "Predict time (s)"),
               ("process_max_mem", "Process peak memory (MB)"), ("predict_max_mem", "Predict peak memory (MB)"),
               ("process_avg_cpu", "Process avg CPU (%)"), ("predict_avg_cpu", "Predict avg CPU (%)")]
APP_METRICS = [("elapsed_time_sec", "Elapsed time (s)"), ("server_max_mem_mb", "Server peak memory (MB)"),
               ("server_avg_cpu_percent", "Server avg CPU (%)"),
               ("browser_max_mem_mb", "Browser peak memory (MB)"), ("browser_avg_cpu_percent", "Browser avg CPU (%)")]
PALETTE = ["#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756"]


def load(results: Path, prefix: str):
    """Return {label: (dataframe, meta)} ordered oldest to newest commit."""
    found = []
    for csv in results.glob(f"{prefix}_*.csv"):
        label = csv.stem[len(prefix) + 1:]
        meta_path = csv.with_suffix(".meta.json")
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        ctime = subprocess.run(["git", "show", "-s", "--format=%ct", meta.get("commit", "HEAD")],
                               cwd=bc.REPO, capture_output=True, text=True).stdout.strip() or "0"
        found.append(((int(ctime), meta.get("dirty", False)), label, pd.read_csv(csv), meta))
    found.sort(key=lambda t: t[0])
    return {label: (df, meta) for _, label, df, meta in found}


def display(label, meta):
    """Short legend name; flags uncommitted working-tree runs."""
    return f"{meta.get('git_describe', label)} (uncommitted)" if meta.get("dirty") else label


def usable_metrics(metrics, df, meta, kind):
    """Drop columns that are missing, or CPU columns from older app runs that summed cumulative CPU."""
    out = []
    for col, name in metrics:
        if col not in df:
            continue
        if kind == "app" and "cpu" in col and not meta.get("params", {}).get("cpu_delta"):
            continue
        out.append((col, name))
    return out


def pm(s, digits=1):
    return f"{s.mean():.{digits}f} ± {s.std():.{digits}f}"


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def html_table(df: pd.DataFrame, num_from=2) -> str:
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    body = "".join("<tr>" + "".join(
        f"<td class={'n' if i >= num_from else 't'}>{html.escape(str(v))}</td>" for i, v in enumerate(row)) + "</tr>"
        for row in df.itertuples(index=False))
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def summary_table(df, group, enzyme, metrics, reads):
    g = df.groupby(group)
    out = pd.DataFrame({"Experiment": list(g.groups), "Enzyme": [x[enzyme].iloc[0] for _, x in g],
                        "Runs": [len(x) for _, x in g]})
    for col, label in metrics:
        out[label] = [pm(x[col], 1 if "mem" not in col else 0) for _, x in g]
    if reads is not None:
        r = reads.groupby("Experiment")[["Control_reads", "Exp_reads"]].mean().astype(int)
        out["Control read pairs"] = [f"{r.loc[int(e), 'Control_reads']:,}" for e in out["Experiment"]]
        out["Experimental read pairs"] = [f"{r.loc[int(e), 'Exp_reads']:,}" for e in out["Experiment"]]
    return out


def figure(data, group, metrics, title, path, ncols=2):
    """One panel per metric; grouped bars (version) per experiment, SD error bars, enzyme tick labels."""
    n = len(metrics)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 3.6 * nrows), squeeze=False)
    labels = list(data)
    first = next(iter(data.values()))[0]
    enz_col = "Enzyme" if "Enzyme" in first else "enzyme"
    exps = sorted(set(first[group].astype(int)))
    names = {int(e): f"{e}\n{first.loc[first[group] == e, enz_col].iloc[0]}" for e in exps}
    width = 0.8 / len(labels)
    for ax, (col, label) in zip(axes.flat, metrics):
        for j, lab in enumerate(labels):
            df, meta = data[lab]
            stats = df.groupby(group)[col].agg(["mean", "std"]).reindex(exps)
            ax.bar([i + j * width for i in range(len(exps))], stats["mean"], width, yerr=stats["std"],
                   capsize=2, color=PALETTE[j % len(PALETTE)], label=display(lab, meta), error_kw={"lw": 0.8})
        ax.set_xticks([i + width * (len(labels) - 1) / 2 for i in range(len(exps))])
        ax.set_xticklabels([names[e] for e in exps], fontsize=8)
        ax.set_title(label, fontsize=11, loc="left", fontweight="bold")
        ax.grid(axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.suptitle(title, fontsize=13, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def comparison_table(cli):
    """Total CLI wall time (process + predict) per experiment for every version, and change vs the first."""
    labels = list(cli)
    rows = None
    for j, lab in enumerate(labels):
        df, meta = cli[lab]
        d = df.assign(total=df.process_time + df.predict_time).groupby("Experiment")
        t = pd.DataFrame({"Experiment": list(d.groups), "Enzyme": [x.Enzyme.iloc[0] for _, x in d],
                          f"{display(lab, meta)} total (s)": [pm(x.total) for _, x in d],
                          f"_m{j}": [x.total.mean() for _, x in d]})
        rows = t if rows is None else rows.merge(t.drop(columns="Enzyme"), on="Experiment")
    if len(labels) > 1:
        rows[f"Change vs {display(labels[0], cli[labels[0]][1])}"] = [
            f"{(b / a - 1) * 100:+.0f}%" for a, b in zip(rows["_m0"], rows[f"_m{len(labels) - 1}"])]
    return rows.drop(columns=[c for c in rows if c.startswith("_m")])


def section(kind, data, results, metrics, group, enzyme, title, png, fig_metrics, reads=None):
    if not data:
        return f"_No {kind} results found in `results/`. Run the benchmark first._\n", {}
    parts, tables = [], {}
    for label, (df, meta) in data.items():
        m = usable_metrics(metrics, df, meta, "app" if kind == "web app" else "cli")
        t = summary_table(df, group, enzyme, m, reads if kind == "web app" else df if "Control_reads" in df else None)
        tables[label] = t
        commit = meta.get("commit", "")[:9]
        dirty = " with uncommitted changes" if meta.get("dirty") else ""
        parts.append(f"#### {display(label, meta)}\n\nCommit `{commit}`{dirty}, run {meta.get('timestamp_utc', '')[:10]}; "
                     f"mean ± SD.\n\n{md_table(t)}\n")
    first_meta = next(iter(data.values()))[1]
    fm = usable_metrics(fig_metrics, next(iter(data.values()))[0], first_meta, "app" if kind == "web app" else "cli")
    figure(data, group, fm, title, results / png)
    parts.append(f"![{title}](results/{png})\n")
    return "\n".join(parts), tables


def system_table(results: Path) -> pd.DataFrame:
    metas = [m for prefix in ("cli", "app") for _, m in load(results, prefix).values() if m]
    if not metas:
        return pd.DataFrame(columns=["Attribute", "Value"])
    s = metas[-1].get("system", {})
    names = {"Model Name": "Model", "Chip": "Chip", "Total Number of Cores": "Cores", "Memory": "Memory",
             "macos": "macOS", "chrome": "Browser", "bbmap": "BBMap", "python": "Python"}
    return pd.DataFrame([(v, s[k]) for k, v in names.items() if k in s], columns=["Attribute", "Value"])


def b64(path):
    return base64.b64encode(Path(path).read_bytes()).decode()


CSS = """
:root{--bg:#fff;--fg:#1f2430;--mut:#667085;--line:#e4e7ec;--head:#f4f6f9;--acc:#4C78A8}
@media(prefers-color-scheme:dark){:root{--bg:#14171c;--fg:#e8eaed;--mut:#9aa3b2;--line:#2c323c;--head:#1c2027;--acc:#7fa6cf}}
body{background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1180px;margin:0 auto;padding:32px 16px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:20px;margin:40px 0 8px;border-bottom:2px solid var(--acc);padding-bottom:4px}
h3{font-size:16px;margin:24px 0 6px}p.s{color:var(--mut);margin:0 0 12px}
.w{overflow-x:auto}table{border-collapse:collapse;font-size:13.5px;margin:8px 0 4px;min-width:60%}
th{background:var(--head);text-align:left;padding:7px 10px;border-bottom:2px solid var(--line);white-space:nowrap}
td{padding:6px 10px;border-bottom:1px solid var(--line);white-space:nowrap}td.n{text-align:right;font-variant-numeric:tabular-nums}
tbody tr:nth-child(even){background:var(--head)}img{max-width:100%;height:auto;margin:10px 0;background:#fff;border-radius:6px}
.note{background:var(--head);border-left:4px solid var(--acc);padding:8px 12px;margin:12px 0;color:var(--mut)}
"""


def build_html(results, cli, app, cli_tables, app_tables, system):
    h = [f"<!doctype html><html lang=en><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
         f"<title>TamiPami performance benchmark</title><style>{CSS}</style><body>",
         "<h1>TamiPami performance benchmark</h1><p class=s>Mean ± SD over repeated runs of each experiment.</p>"]
    dirty = [display(lab, m) for lab, (_, m) in {**cli, **app}.items() if m.get("dirty")]
    if dirty:
        h.append("<div class=note>Runs marked “uncommitted” used the working tree, which differs from the named commit. "
                 "Commit and tag before publishing numbers.</div>")
    if cli:
        h.append("<h2>Command line</h2><h3>Total time per dataset (process + predict)</h3><div class=w>"
                 + html_table(comparison_table(cli)) + "</div>")
        h.append(f"<img alt='CLI benchmark' src='data:image/png;base64,{b64(results / 'cli_benchmark.png')}'>")
        for lab, t in cli_tables.items():
            h.append(f"<h3>{html.escape(display(lab, cli[lab][1]))}</h3><div class=w>{html_table(t)}</div>")
    if app:
        h.append("<h2>Web application</h2>")
        h.append(f"<img alt='App benchmark' src='data:image/png;base64,{b64(results / 'app_benchmark.png')}'>")
        for lab, t in app_tables.items():
            h.append(f"<h3>{html.escape(display(lab, app[lab][1]))}</h3><div class=w>{html_table(t)}</div>")
        if not all(m.get("params", {}).get("cpu_delta") for _, m in app.values()):
            h.append("<div class=note>CPU columns are omitted for app runs recorded before the sampler measured CPU "
                     "<em>during</em> each submission; the earlier values summed the long-lived server's cumulative CPU time. "
                     "Server memory is the peak resident size of the one server process, which accumulates across the "
                     "runs in a session.</div>")
    h.append("<h2>Test system</h2><div class=w>" + html_table(system, 1) + "</div></body></html>")
    return "".join(h)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=str(bc.RESULTS_DIR))
    ap.add_argument("--readme", default=str(bc.BENCH_DIR / "readme.md"))
    args = ap.parse_args()
    results = Path(args.results)

    cli = load(results, "cli")
    app = load(results, "app")
    cli_md, cli_t = section("CLI", cli, results, CLI_METRICS, "Experiment", "Enzyme", "Command line benchmark",
                            "cli_benchmark.png", CLI_METRICS[:4])
    reads = list(cli.values())[-1][0] if cli else None  # app runs reuse the CLI read counts (same inputs)
    app_md, app_t = section("web app", app, results, APP_METRICS, "experiment", "enzyme", "Web application benchmark",
                            "app_benchmark.png", APP_METRICS, reads)
    system = system_table(results)
    sys_md = md_table(system) if len(system) else ""
    (results / "cli_summary.md").write_text(cli_md)
    (results / "app_summary.md").write_text(app_md)
    (results / "system.md").write_text(sys_md)
    (results / "report.html").write_text(build_html(results, cli, app, cli_t, app_t, system))

    block = (f"## Test system configuration\n\nRetrieved automatically when the benchmarks ran.\n\n{sys_md}\n\n"
             f"## Results\n\nGenerated by `summarize_benchmarks.py`; a styled version is in `results/report.html`.\n\n"
             f"### CLI results\n\n{cli_md}\n### Web application results\n\n{app_md}")
    readme = Path(args.readme)
    new = re.sub(r"(<!-- BEGIN AUTO -->\n).*?(<!-- END AUTO -->)",
                 lambda m: m.group(1) + block + "\n" + m.group(2), readme.read_text(), flags=re.S)
    readme.write_text(new)
    print(f"Updated {readme} and {results}")


if __name__ == "__main__":
    main()
