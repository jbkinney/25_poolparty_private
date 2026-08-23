"""Builds figure.drawio for the VaLiAnT comparison.

Cell vocabulary from ../../mechanism_figure/fig_s2.drawio; width budget from
../README.md.  Edit this script, not the XML.

Panel A is drawn in GENOMIC order, left to right, which is the order r1/r2/r3
are defined in: valiant/loaders/targeton_config.py assigns r1 = region_2
.get_before(n) and r3 = .get_after(n), and uint_range.py's get_before/get_after
are pure coordinate arithmetic with no strand term.  So on this minus-strand
targeton the 1-nt orphan carrying phase = 1 sits at the RIGHT edge of r2
(43,115,779, the transcript 5' end) and the 2-nt orphan at the LEFT
(43,115,726-727).  An earlier ASCII layout had these two annotations swapped.
"""
import os, re
from xml.sax.saxutils import escape

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
SEG  = "text;html=1;strokeColor=none;fillColor=none;align=center;verticalAlign=middle;fontSize=7;"
SM   = "text;html=1;strokeColor=none;fillColor=none;align=center;verticalAlign=middle;fontSize=6.5;fontColor=#666666;"

cells, _n = [], [0]
def cell(style, value, x, y, w, h, tag=""):
    _n[0] += 1
    cells.append((tag, f'        <mxCell id="c{_n[0]}" parent="1" style="{style}" value="{value}" vertex="1">\n'
                       f'          <mxGeometry x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" as="geometry" />\n'
                       f'        </mxCell>', (x, y, w, h)))

def code_html(lines, colour_strings=True):
    out = []
    for ln in lines:
        if not ln.strip():
            out.append("<div><br></div>"); continue
        if colour_strings:
            parts, last = [], 0
            for m in re.finditer(r'"[^"]*"', ln):
                parts.append(escape(ln[last:m.start()]))
                parts.append(f'<span style="color: {STR_COL};">{escape(m.group(0))}</span>')
                last = m.end()
            parts.append(escape(ln[last:]))
            body = "".join(parts)
        else:
            body = escape(ln)
        lead = len(ln) - len(ln.lstrip(" "))
        if lead:
            body = "&nbsp;" * lead + body.lstrip(" ")
        out.append(f'<div><font style="font-size: 8.3px;">{body}</font></div>')
    return "".join(out)

W = 670.0
# ---------------------------------------------------------------- panel A
BX, BW, BH = 120.0, 430.0, 16.0
TOT = 104                                   # r1 25 + r2 54 + r3 25
def nt(n): return BX + n / TOT * BW
R1, R2A, R2B, R3 = 0, 25, 79, 104           # genomic order, left to right
O_L, O_R = 27, 78                           # r2 interior: 2 nt | 51 nt | 1 nt

Y_T, Y_COORD, Y_R2, Y_BAR = 0.0, 16.0, 32.0, 46.0
Y_ENDS, Y_ANN, Y_S1, Y_S2 = Y_BAR + BH + 2, Y_BAR + BH + 18, 108.0, 122.0

cell(LBL, "A", 0, Y_T, 14, 14)
cell(TXT, "One library: BRCA1 exon 2, three mutation types, P5/P7 adaptors", 18, Y_T, 340, 14)
cell(TXT, "chr17:43,115,634-43,115,878   minus strand, 245 nt", BX, Y_COORD, 300, 12)
cell(CTR, "r2 = exon 2, 54 nt", nt(R2A), Y_R2, nt(R2B) - nt(R2A), 12)

cell("rounded=0;whiteSpace=wrap;html=1;fillColor=#F7F7F7;strokeColor=#999999;", "",
     nt(R1), Y_BAR, nt(R2A) - nt(R1), BH, "bar")
cell("rounded=0;whiteSpace=wrap;html=1;fillColor=#F7F7F7;strokeColor=#999999;", "",
     nt(R2B), Y_BAR, nt(R3) - nt(R2B), BH, "bar")
cell("rounded=0;whiteSpace=wrap;html=1;fillColor=#DAE8FC;strokeColor=#6C8EBF;", "",
     nt(R2A), Y_BAR, nt(R2B) - nt(R2A), BH, "bar")
for a, b in ((R2A, O_L), (O_R, R2B)):        # the two orphaned-base blocks
    cell("rounded=0;whiteSpace=wrap;html=1;fillColor=#F8CECC;strokeColor=#B85450;", "",
         nt(a), Y_BAR, nt(b) - nt(a), BH, "bar")
cell(SEG, "r1  25 nt", nt(R1), Y_BAR, nt(R2A) - nt(R1), BH, "bar")
cell(SEG, "r3  25 nt", nt(R2B), Y_BAR, nt(R3) - nt(R2B), BH, "bar")
cell(SEG, "17 codons, 51 nt", nt(O_L), Y_BAR, nt(O_R) - nt(O_L), BH, "bar")

cell(SM, "transcript 3&#39; end", nt(R1) - 4, Y_ENDS, 76, 10)
cell(SM, "transcript 5&#39; end", nt(R3) - 72, Y_ENDS, 76, 10)
for x, txt, wid in ((nt(R2A) + (nt(O_L) - nt(R2A)) / 2,
                     "2 nt &#8211; begins a codon completed in exon 3", 190),
                    (nt(O_R) + (nt(R2B) - nt(O_R)) / 2,
                     "1 nt &#8211; completes a codon begun in exon 1;  phase = 1", 236)):
    cell("endArrow=none;html=1;strokeColor=#B85450;strokeWidth=0.5;", "",
         x - 0.5, Y_BAR + BH, 1, Y_ANN - Y_BAR - BH, "leader")
    cell(CTR + "fontColor=#B85450;", txt, x - wid / 2, Y_ANN, wid, 12)

cell(TXT, "17 codons &#215; (19 aa + 1 in-frame deletion + 1 stop)  =  357 oligos",
     BX, Y_S1, 360, 12)
cell(TXT, "P5 20 nt  +  245 nt  +  P7 21 nt  =  286 nt  (283 where a codon is deleted)", BX, Y_S2, 400, 12)

# ---------------------------------------------------------------- panels B, C
B = ['targeton.tsv',
     '   ref_chr        chr17',
     '   ref_strand     -',
     '   ref_start      43115634',
     '   ref_end        43115878',
     '   r2_start       43115726',
     '   r2_end         43115779',
     '   ext_vector     "25,25"',
     '   action_vector  "(),(aa,inframe,stop),()"',
     '',
     'ENST00000357654.9.gtf',
     '   CDS  43115726  43115779  .  -  1',
     '',
     'chr17.fa   81 MB',
     '',
     '$ valiant sge targeton.tsv ../ref/chr17.fa out \\',
     "      'homo sapiens' 'GRCh38' \\",
     '      --revcomp-minus-strand \\',
     '      --adaptor-5 ... --adaptor-3 ... \\',
     '      --gff ENST00000357654.9.gtf',
     '',
     '3 input files + one command line']
C = ['ADAPTOR_5 = "AATGATACGGCGACCACCGA"',
     'ADAPTOR_3 = "TCGTATGCCGTCTTCTGCTTG"',
     '',
     'orf = pp.from_fasta(',
     '    "chr17.fa", ("chr17", 43115633, 43115878, "-")',
     ').annotate_orf(region_name="cds", extent=(100, 151), frame=1)',
     '',
     'aa = orf.mutagenize_orf(',
     '    region="cds", num_mutations=1,',
     '    mutation_type="missense_only_first", mode="sequential")',
     '',
     'inframe = orf.deletion_scan_orf(',
     '    deletion_codons=1, deletion_marker=None,',
     '    region="cds", mode="sequential")',
     '',
     'stop = orf.insertion_scan_orf(',
     '    "TGA", region="cds", replace=True, mode="sequential")',
     '',
     'library = pp.join([ADAPTOR_5, pp.stack([aa, inframe, stop]),',
     '                   ADAPTOR_3])']

BW_, CW_ = 56 * CHAR, 70 * CHAR
XB, XC = 0.0, W - CW_
Y_PANEL = 146.0
Y_CODE = Y_PANEL + 18

for letter, name, code, x, wid, colour in (("B", "VaLiAnT 4.0.0", B, XB, BW_, False),
                                           ("C", "PoolParty",     C, XC, CW_, True)):
    cell(LBL, letter, x, Y_PANEL, 14, 14)
    cell(NAME, name, x + 18, Y_PANEL, 160, 14)
    cell(CODE, escape(code_html(code, colour), {'"': "&quot;"}),
         x, Y_CODE, wid, len(code) * LINE + 12, "code")

# phase -> frame=1, the one value transcribed by hand
def at(col, idx, x0): return (x0 + 8 + col * CHAR, Y_CODE + 6 + idx * LINE + LINE / 2)
px, py = at(len(B[11]) - 1, 11, XB)
fx, fy = at(C[5].index("frame=1"), 5, XC)
_n[0] += 1
cells.append(("edge",
    f'        <mxCell id="c{_n[0]}" parent="1" edge="1" '
    f'style="edgeStyle=orthogonalEdgeStyle;html=1;rounded=1;strokeColor=#B85450;'
    f'strokeWidth=0.75;endArrow=blockThin;endFill=1;endSize=4;exitX=0.5;exitY=1;'
    f'fontSize=7;fontColor=#B85450;labelBackgroundColor=#FFFFFF;" '
    f'value="transcribed by hand">\n'
    f'          <mxGeometry relative="1" as="geometry">\n'
    f'            <mxPoint x="{px:.1f}" y="{py:.1f}" as="sourcePoint" />\n'
    f'            <mxPoint x="{fx:.1f}" y="{fy:.1f}" as="targetPoint" />\n'
    f'            <Array as="points"><mxPoint x="{px:.1f}" y="{fy - 26:.1f}" />'
    f'<mxPoint x="{fx:.1f}" y="{fy - 26:.1f}" /></Array>\n'
    f'          </mxGeometry>\n        </mxCell>', None))

H = Y_CODE + max(len(B), len(C)) * LINE + 12
xml = ('<mxfile host="app.diagrams.net">\n'
       '  <diagram name="Page-1" id="valiant-comparison">\n'
       f'    <mxGraphModel dx="900" dy="700" grid="1" gridSize="10" guides="1" tooltips="1" '
       f'connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="{W:.0f}" '
       f'pageHeight="{H:.0f}" math="0" shadow="0">\n      <root>\n'
       '        <mxCell id="0" />\n        <mxCell id="1" parent="0" />\n'
       + "\n".join(c[1] for c in cells) + '\n      </root>\n    </mxGraphModel>\n  </diagram>\n</mxfile>\n')
open(DEST, "w", encoding="utf-8").write(xml)

# ---- checks -----------------------------------------------------------------
def hit(a, b):
    return a[0] < b[0]+b[2] and b[0] < a[0]+a[2] and a[1] < b[1]+b[3] and b[1] < a[1]+a[3]
boxed = [(t, g) for t, _, g in cells if g and t not in ("bar", "leader")]
clash = [(boxed[i][0] or "?", boxed[j][0] or "?")
         for i in range(len(boxed)) for j in range(i+1, len(boxed))
         if hit(boxed[i][1], boxed[j][1])]
print(f"figure.drawio: {len(cells)} cells, page {W:.0f} x {H:.0f}")
print(f"  aspect {W/H:.2f}:1  ->  at 482 pt wide, {482*H/W:.0f} pt / {482*H/W/72*25.4:.0f} mm tall")
print(f"  panel B {BW_/CHAR:.0f} chars, panel C {CW_/CHAR:.0f} chars, full width {W/CHAR:.0f}")
print(f"  longest line: B {max(len(l) for l in B)}, C {max(len(l) for l in C)}")
print(f"  overlaps (excl. bar/leaders): {clash if clash else 'none'}")


# ---- arithmetic, asserted so the figure cannot state a wrong total ----------
N_CODONS, N_AA, N_DEL, N_STOP = 17, 19, 1, 1
assert (O_L - R2A) + (O_R - O_L) + (R2B - O_R) == 54          # 2 + 51 + 1
assert (O_R - O_L) == N_CODONS * 3                            # 51 nt = 17 codons
assert (R2A - R1) + 54 + (R3 - R2B) == TOT == 104
assert N_CODONS * (N_AA + N_DEL + N_STOP) == 357
assert 20 + 245 + 21 == 286 and 286 - 3 == 283

# ---- emit the ASCII record in figure.md from the same numbers ---------------
WCH, BCH, CCH, GUT = 130, 56, 70, 4

def box(title, content, w, indent=2):
    pad = " " * indent
    return (["\u250c" + "\u2500" * (w - 2) + "\u2510",
             f"\u2502 {title}".ljust(w - 1) + "\u2502",
             "\u2502".ljust(w - 1) + "\u2502"]
            + [("\u2502" + pad + c).ljust(w - 1) + "\u2502" for c in content]
            + ["\u2514" + "\u2500" * (w - 2) + "\u2518"])

A_ASCII = [
    "chr17:43,115,634-43,115,878   minus strand, 245 nt",
    "\u250c" + "\u2500" * 10 + "\u252c" + "\u2500" * 3 + "\u252c" + "\u2500" * 30
        + "\u252c" + "\u2500" * 3 + "\u252c" + "\u2500" * 10 + "\u2510",
    "\u2502 r1 25 nt \u2502 2 \u2502      17 codons, 51 nt        \u2502 1 \u2502 r3 25 nt \u2502",
    "\u2514" + "\u2500" * 10 + "\u2534" + "\u2500" * 3 + "\u2534" + "\u2500" * 30
        + "\u2534" + "\u2500" * 3 + "\u2534" + "\u2500" * 10 + "\u2518",
    " transcript 3'  \u25b2" + " " * 30 + "\u25b2   transcript 5'",
    " " * 16 + "\u2502" + " " * 30 + "\u2502",
    " " * 6 + "begins a codon completed         completes a codon begun in",
    " " * 11 + "in exon 3                        exon 1;  phase = 1",
    "",
    "17 codons \u00d7 (19 aa + 1 in-frame deletion + 1 stop)  =  357 oligos",
    "P5 20 nt + 245 nt + P7 21 nt = 286 nt  (283 where a codon is deleted)",
]

out = box("(A)  One library: BRCA1 exon 2, three mutation types, P5/P7 adaptors",
          A_ASCII, WCH, indent=4)
b, c = box("(B)  VaLiAnT 4.0.0", B, BCH), box("(C)  PoolParty", C, CCH)
def pad_to(lines, n, w):
    """Insert filler rows before the closing border, not after it."""
    if len(lines) >= n:
        return lines
    return lines[:-1] + ["\u2502".ljust(w - 1) + "\u2502"] * (n - len(lines)) + [lines[-1]]

n = max(len(b), len(c))
b, c = pad_to(b, n, BCH), pad_to(c, n, CCH)
out += [lb + " " * GUT + lc for lb, lc in zip(b, c)]
block = "\n".join(out)

md_path = os.path.join(os.path.dirname(DEST), "figure.md")
md = open(md_path, encoding="utf-8").read()
i = md.index("## Layout\n\n```\n") + len("## Layout\n\n```\n")
j = md.index("\n```\n\n---\n\n## Panel notes")
open(md_path, "w", encoding="utf-8").write(md[:i] + block + md[j:])
print(f"  figure.md layout regenerated: {len(out)} lines, widths {sorted({len(l) for l in out})}")
