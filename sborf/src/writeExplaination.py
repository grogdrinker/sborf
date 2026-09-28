import math
import tempfile
import webbrowser
import os
import sys


def _display_available():
    # A browser can realistically be opened on macOS/Windows, or on a
    # Linux/Unix session that exposes an X11 / Wayland display.
    if sys.platform == "darwin" or sys.platform.startswith("win"):
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def print_codons_terminal(sequence, scores, max_abs=1.0):
    # Prints the per-codon scores as an aligned table on stdout.
    # Uses ANSI background colors only when writing to an interactive terminal,
    # mirroring the green(positive)/red(negative) scheme of the HTML output.
    # Scores are expected in [-1, 1] (fixed color scale, max_abs = 1).
    use_color = sys.stdout.isatty()

    def colorize(text, score):
        if not use_color or score == 0:
            return text
        intensity = min(abs(score) / max_abs, 1)
        if score > 0:
            code = 42 if intensity > 0.5 else 102  # green background
        else:
            code = 41 if intensity > 0.5 else 101  # red background
        return "\033[{}m{}\033[0m".format(code, text)

    print("{:>4}  {:<6} {:>10}".format("Pos", "Codon", "Score"))
    print("-" * 24)
    for i, (codon, score) in enumerate(zip(sequence, scores)):
        line = "{:>4}  {:<6} {:>10.4f}".format(i, codon, score)
        print(colorize(line, score))


def visualize_codons(sequence, scores, max_abs=1.0):
    # sequence: list of codons, e.g. ["AAA","AAC","AAA","UGG",...]
    # scores:   list of floats of same length, in [-1, 1] (fixed color scale, max_abs = 1)

    def score_to_color(score):
        x = abs(score) / max_abs
        x = min(max(x, 0), 1)
        if score > 0:
            r = 255 - int(120 * x)
            g = 255
            b = 255 - int(120 * x)
        elif score < 0:
            r = 255
            g = 255 - int(120 * x)
            b = 255 - int(120 * x)
        else:
            return "rgb(255,255,255)"
        return f"rgb({r},{g},{b})"

    html = [
        "<html><body>",
        "<h2>Codon Score Visualization</h2>",
        "<table border='1' cellpadding='6' style='border-collapse:collapse'>",
        "<tr><th>Position</th><th>Codon</th><th>Score</th></tr>"
    ]

    for i, (codon, score) in enumerate(zip(sequence, scores)):
        color = score_to_color(score)
        html.append(
            f"<tr>"
            f"<td style='background:{color}'>{i}</td>"
            f"<td style='background:{color}'>{codon}</td>"
            f"<td style='background:{color}'>{score}</td>"
            f"</tr>"
        )

    html.extend(["</table>", "</body></html>"])
    document = "\n".join(html)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as f:
        f.write(document.encode("utf-8"))
        path = os.path.abspath(f.name)

    if _display_available():
        webbrowser.open("file://" + path)
        print("HTML visualization opened in your browser (also saved to: " + path + ")")
    else:
        # Headless session: don't silently try to launch a browser.
        # Show the scores in the terminal and point to the saved HTML file.
        print_codons_terminal(sequence, scores, max_abs)
        print("HTML visualization saved to: " + path)


def write_codons_visualization(sequence, scores,outfile, max_abs=1.0):
    # sequence: list of codons, e.g. ["AAA","AAC","AAA","UGG",...]
    # scores:   list of floats of same length, in [-1, 1] (fixed color scale, max_abs = 1)

    def score_to_color(score):
        x = abs(score) / max_abs
        x = min(max(x, 0), 1)
        if score > 0:
            r = 255 - int(120 * x)
            g = 255
            b = 255 - int(120 * x)
        elif score < 0:
            r = 255
            g = 255 - int(120 * x)
            b = 255 - int(120 * x)
        else:
            return "rgb(255,255,255)"
        return f"rgb({r},{g},{b})"

    html = [
        "<html><body>",
        "<h2>Codon Score Visualization</h2>",
        "<table border='1' cellpadding='6' style='border-collapse:collapse'>",
        "<tr><th>Position</th><th>Codon</th><th>Score</th></tr>"
    ]

    for i, (codon, score) in enumerate(zip(sequence, scores)):
        color = score_to_color(score)
        html.append(
            f"<tr>"
            f"<td style='background:{color}'>{i}</td>"
            f"<td style='background:{color}'>{codon}</td>"
            f"<td style='background:{color}'>{score}</td>"
            f"</tr>"
        )

    html.extend(["</table>", "</body></html>"])
    document = "\n".join(html)

    with open(outfile, "w", encoding="utf-8") as f:
        f.write(document)