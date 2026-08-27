# Builds figure.drawio for the tangermeme comparison, using the cell vocabulary of
# ../../mechanism_figure/fig_s2.drawio.  Spec: code_comparison/README.md.
import re
from xml.sax.saxutils import escape

import os
DEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figure.drawio")
CHAR, LINE = 5.16, 11.0
STR_COL = "rgb(163, 21, 21)"

CODE = ("text;html=1;strokeColor=none;fillColor=none;align=left;verticalAlign=top;"
        "whiteSpace=wrap;rounded=0;fontSize=9;fontFamily=Consolas;spacing=5;"
        "spacingLeft=8;spacingTop=6;")
LBL  = "text;html=1;strokeColor=none;fillColor=none;align=left;verticalAlign=middle;fontSize=10;fontStyle=1;"
NAME = "text;html=1;strokeColor=none;fillColor=none;align=left;verticalAlign=middle;fontSize=9;"
TXT  = "text;html=1;strokeColor=none;fillColor=none;align=left;verticalAlign=middle;fontSize=8;"
CTR  = "text;html=1;strokeColor=none;fillColor=none;align=center;verticalAlign=middle;fontSize=8;"
TICK = "text;html=1;strokeColor=none;fillColor=none;align=center;verticalAlign=middle;fontSize=6.5;"

cells, _n = [], [0]
def cell(style, value, x, y, w, h):
    _n[0] += 1
    cells.append(f'        <mxCell id="c{_n[0]}" parent="1" style="{style}" value="{value}" vertex="1">\n'
                 f'          <mxGeometry x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" as="geometry" />\n'
                 f'        </mxCell>')

def code_html(lines):
    out = []
    for ln in lines:
        if not ln.strip():
            out.append("<div><br></div>"); continue
        parts, last = [], 0
        for m in re.finditer(r'"[^"]*"', ln):          # only strings are coloured
            parts.append(escape(ln[last:m.start()]))
            parts.append(f'<span style="color: {STR_COL};">{escape(m.group(0))}</span>')
            last = m.end()
        parts.append(escape(ln[last:]))
        body = "".join(parts)
        lead = len(ln) - len(ln.lstrip(" "))
        if lead:
            body = "&nbsp;" * lead + body.lstrip(" ")
        out.append(f'<div><font style="font-size: 8.3px;">{body}</font></div>')
    return "".join(out)

W = 670.0
BX, BW, BH = 135.0, 400.0, 15.0
def bp(p): return BX + p / 201.0 * BW

# ---- panel A -------------------------------------------------- vertical bands
Y_DESC, Y_NOTE, Y_ARR, Y_BAR = 0.0, 22.0, 34.0, 48.0
Y_TICK, Y_SUM = Y_BAR + BH + 1, Y_BAR + BH + 19

cell(LBL, "A", 0, Y_DESC, 14, 14)
cell(TXT, "201 bp region containing the canonical 5&#39; splice site of SMN2 exon 7",
     18, Y_DESC, 400, 14)
cell(CTR + "fontColor=#A50040;", "canonical 5&#39; splice site, not scanned",
     bp(100) - 90, Y_NOTE, 180, 12)
cell("endArrow=none;html=1;strokeColor=#A50040;strokeWidth=1;", "",
     bp(100) - 1, Y_ARR, 2, Y_BAR - Y_ARR)
cell("rounded=0;whiteSpace=wrap;html=1;fillColor=#F7F7F7;strokeColor=#999999;", "",
     BX, Y_BAR, BW, BH)
for a, b in ((51, 89), (107, 167)):
    cell("rounded=0;whiteSpace=wrap;html=1;fillColor=#DAE8FC;strokeColor=#6C8EBF;", "",
         bp(a), Y_BAR, bp(b) - bp(a), BH)
cell("rounded=0;whiteSpace=wrap;html=1;fillColor=#A50040;strokeColor=none;", "",
     bp(99), Y_BAR, bp(101) - bp(99), BH)
for p in (1, 51, 89, 107, 167, 201):
    cell(TICK, str(p), bp(p) - 12, Y_TICK, 24, 9)
cell(TXT, "100 scan positions", BX + BW + 8, Y_BAR, 120, BH)
cell(TXT, "2,000 cryptic 5&#39; splice sites  &#215;  100 positions  =  200,000 sequences",
     BX, Y_SUM, 400, 14)

# ---- panels B and C ---------------------------------------------------------
B = ['cryptic_pool = pp.from_seqs(',
     '    cryptic_sites, mode="sequential",',
     '    cards={"seq": "cryptic_sequence"})',
     'library = pp.from_seq(target_region).replacement_scan(',
     '    cryptic_pool, positions=SCAN_POSITIONS,',
     '    mode="sequential",',
     '    cards={"start": "cryptic_position"})',
     'df = library.generate_library()']
C = ['target_ohe = (one_hot_encode(target_region)',
     '              .unsqueeze(0)',
     '              .expand(len(cryptic_sites), -1, -1))',
     'cryptic_ohe = torch.stack(',
     '    [one_hot_encode(site) for site in cryptic_sites])',
     'library_ohe = torch.stack(',
     '    [substitute(target_ohe, cryptic_ohe, start=p)',
     '     for p in SCAN_POSITIONS])']

Y_PANEL, CW = Y_SUM + 30, 330.0
for letter, name, code, x in (("B", "PoolParty", B, 0.0), ("C", "tangermeme", C, W - CW)):
    cell(LBL, letter, x, Y_PANEL, 14, 14)
    cell(NAME, name, x + 18, Y_PANEL, 140, 14)
    cell(CODE, escape(code_html(code), {'"': "&quot;"}),
         x, Y_PANEL + 18, CW, len(code) * LINE + 12)

H = Y_PANEL + 18 + len(B) * LINE + 12
xml = ('<mxfile host="app.diagrams.net">\n'
       '  <diagram name="Page-1" id="tangermeme-comparison">\n'
       f'    <mxGraphModel dx="900" dy="600" grid="1" gridSize="10" guides="1" tooltips="1" '
       f'connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="{W:.0f}" '
       f'pageHeight="{H:.0f}" math="0" shadow="0">\n      <root>\n'
       '        <mxCell id="0" />\n        <mxCell id="1" parent="0" />\n'
       + "\n".join(cells) + '\n      </root>\n    </mxGraphModel>\n  </diagram>\n</mxfile>\n')
open(DEST, "w", encoding="utf-8").write(xml)

# ---- overlap check ----------------------------------------------------------
boxes = [(float(m[0]), float(m[1]), float(m[2]), float(m[3]))
         for m in re.findall(r'x="([-\d.]+)" y="([-\d.]+)" width="([\d.]+)" height="([\d.]+)"', xml)]
def hit(a, b):
    return (a[0] < b[0]+b[2] and b[0] < a[0]+a[2] and a[1] < b[1]+b[3] and b[1] < a[1]+a[3])
# the bar's shaded sub-regions are meant to sit inside the bar, so exempt that band
bar_band = lambda r: abs(r[1] - Y_BAR) < 0.5
clashes = [(i, j) for i in range(len(boxes)) for j in range(i+1, len(boxes))
           if hit(boxes[i], boxes[j]) and not (bar_band(boxes[i]) and bar_band(boxes[j]))]
print(f"figure.drawio: {len(cells)} cells, page {W:.0f} x {H:.0f}")
print(f"  aspect {W/H:.2f}:1  ->  at 482 pt wide, {482*H/W:.0f} pt / {482*H/W/72*25.4:.0f} mm tall")
print(f"  code column {CW:.0f} = {CW/CHAR:.0f} chars;  full width {W/CHAR:.0f} chars")
print(f"  overlaps outside the bar band: {clashes if clashes else 'none'}")
