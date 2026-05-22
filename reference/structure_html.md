# Embed a tRNA structure SVG as centered HTML

Reads an SVG file (typically produced by
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md))
and wraps it in a centering `<div>`, returning an
[`htmltools::HTML()`](https://rstudio.github.io/htmltools/reference/HTML.html)
object suitable for use in R Markdown or Quarto documents.

## Usage

``` r
structure_html(svg_path)
```

## Arguments

- svg_path:

  Path to an SVG file, typically the return value of
  [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Value

An
[`htmltools::HTML()`](https://rstudio.github.io/htmltools/reference/HTML.html)
object.

## Examples

``` r
# \donttest{
svg <- plot_tRNA_structure("tRNA-Glu-TTC", "Escherichia coli")
structure_html(svg)
#> <div style="text-align: center;"><?xml version="1.0" encoding="UTF-8"?>
#> <svg xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.0" id="svg2403" xml:space="preserve" viewBox="0 0 157.727 191.837">
#> <text x="0" y="8.616" id="text1000">
#>   <tspan x="0" y="8.616" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="9" id="tspan1001">tRNA-Glu-TTC</tspan>
#> </text>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 81.6597 49.0485 L 85.5477 49.0485 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 81.6597 56.6085 L 85.5477 56.6085 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 81.6597 64.1685 L 85.5477 64.1685 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 81.6597 71.7285 L 85.5477 71.7285 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 81.6597 79.2885 L 85.5477 79.2885 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 81.6597 86.8485 L 85.5477 86.8485 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 81.6597 94.4085 L 85.5477 94.4085 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 60.0967 110.885 L 60.0967 114.773 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 52.5367 110.885 L 52.5367 114.773 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 44.9767 110.885 L 44.9767 114.773 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 37.4167 110.885 L 37.4167 114.773 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 74.7728 129.42 L 78.6608 129.42 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 74.7728 136.98 L 78.6608 136.98 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 74.7728 144.54 L 78.6608 144.54 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 74.7728 152.1 L 78.6608 152.1 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 74.7728 159.66 L 78.6608 159.66 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 95.0847 107.803 L 95.0847 103.915 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 102.645 107.803 L 102.645 103.915 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 110.205 107.803 L 110.205 103.915 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 117.765 107.803 L 117.765 103.915 "/>
#> <path fill="none" stroke="#999999" stroke-width="1.44" stroke-linecap="butt" stroke-linejoin="miter" d="M 125.325 107.803 L 125.325 103.915 "/>
#> <text x="74.5662" y="51.741" id="text1002">
#>   <tspan x="74.5662" y="51.741" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1003">G</tspan>
#> </text>
#> <text x="74.7762" y="59.301" id="text1004">
#>   <tspan x="74.7762" y="59.301" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1005">U</tspan>
#> </text>
#> <text x="74.7762" y="66.861" id="text1006">
#>   <tspan x="74.7762" y="66.861" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1007">C</tspan>
#> </text>
#> <text x="74.7762" y="74.421" id="text1008">
#>   <tspan x="74.7762" y="74.421" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1009">C</tspan>
#> </text>
#> <text x="74.7762" y="81.981" id="text1010">
#>   <tspan x="74.7762" y="81.981" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1011">C</tspan>
#> </text>
#> <text x="74.7762" y="89.541" id="text1012">
#>   <tspan x="74.7762" y="89.541" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1013">C</tspan>
#> </text>
#> <text x="74.7762" y="97.101" id="text1014">
#>   <tspan x="74.7762" y="97.101" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1015">U</tspan>
#> </text>
#> <text x="67.345" y="98.4906" id="text1016">
#>   <tspan x="67.345" y="98.4906" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1017">U</tspan>
#> </text>
#> <text x="61.1733" y="102.857" id="text1018">
#>   <tspan x="61.1733" y="102.857" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1019">C</tspan>
#> </text>
#> <text x="57.1792" y="109.402" id="text1020">
#>   <tspan x="57.1792" y="109.402" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1021">G</tspan>
#> </text>
#> <text x="49.8292" y="109.402" id="text1022">
#>   <tspan x="49.8292" y="109.402" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1023">U</tspan>
#> </text>
#> <text x="42.2692" y="109.402" id="text1024">
#>   <tspan x="42.2692" y="109.402" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1025">C</tspan>
#> </text>
#> <text x="34.7092" y="109.402" id="text1026">
#>   <tspan x="34.7092" y="109.402" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1027">U</tspan>
#> </text>
#> <text x="29.9567" y="103.695" id="text1028">
#>   <tspan x="29.9567" y="103.695" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1029">A</tspan>
#> </text>
#> <text x="22.3554" y="101.344" id="text1030">
#>   <tspan x="22.3554" y="101.344" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1031">G</tspan>
#> </text>
#> <text x="15.3987" y="103.015" id="text1032">
#>   <tspan x="15.3987" y="103.015" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1033">A</tspan>
#> </text>
#> <text x="9.51341" y="108.235" id="text1034">
#>   <tspan x="9.51341" y="108.235" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1035">G</tspan>
#> </text>
#> <text x="7.5" y="115.522" id="text1036">
#>   <tspan x="7.5" y="115.522" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1037">G</tspan>
#> </text>
#> <text x="9.72341" y="122.809" id="text1038">
#>   <tspan x="9.72341" y="122.809" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1039">C</tspan>
#> </text>
#> <text x="15.1924" y="128.028" id="text1040">
#>   <tspan x="15.1924" y="128.028" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1041">C</tspan>
#> </text>
#> <text x="22.5654" y="129.699" id="text1042">
#>   <tspan x="22.5654" y="129.699" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1043">C</tspan>
#> </text>
#> <text x="29.9567" y="127.348" id="text1044">
#>   <tspan x="29.9567" y="127.348" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1045">A</tspan>
#> </text>
#> <text x="34.4992" y="121.642" id="text1046">
#>   <tspan x="34.4992" y="121.642" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1047">G</tspan>
#> </text>
#> <text x="42.0592" y="121.642" id="text1048">
#>   <tspan x="42.0592" y="121.642" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1049">G</tspan>
#> </text>
#> <text x="50.0354" y="121.642" id="text1050">
#>   <tspan x="50.0354" y="121.642" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1051">A</tspan>
#> </text>
#> <text x="57.3892" y="121.642" id="text1052">
#>   <tspan x="57.3892" y="121.642" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1053">C</tspan>
#> </text>
#> <text x="61.8028" y="127.923" id="text1054">
#>   <tspan x="61.8028" y="127.923" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1055">A</tspan>
#> </text>
#> <text x="67.8893" y="132.113" id="text1056">
#>   <tspan x="67.8893" y="132.113" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1057">C</tspan>
#> </text>
#> <text x="67.8893" y="139.673" id="text1058">
#>   <tspan x="67.8893" y="139.673" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1059">C</tspan>
#> </text>
#> <text x="67.6793" y="147.233" id="text1060">
#>   <tspan x="67.6793" y="147.233" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1061">G</tspan>
#> </text>
#> <text x="67.8893" y="154.793" id="text1062">
#>   <tspan x="67.8893" y="154.793" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1063">C</tspan>
#> </text>
#> <text x="67.8893" y="162.353" id="text1064">
#>   <tspan x="67.8893" y="162.353" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1065">C</tspan>
#> </text>
#> <text x="63.0018" y="168.12" id="text1066">
#>   <tspan x="63.0018" y="168.12" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1067">C</tspan>
#> </text>
#> <text x="62.5989" y="175.67" id="text1068">
#>   <tspan x="62.5989" y="175.67" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1069">U</tspan>
#> </text>
#> <text x="66.8446" y="181.925" id="text1070">
#>   <tspan x="66.8446" y="181.925" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1071">U</tspan>
#> </text>
#> <text x="74.0093" y="184.337" id="text1072">
#>   <tspan x="74.0093" y="184.337" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1073">U</tspan>
#> </text>
#> <text x="81.1739" y="181.925" id="text1074">
#>   <tspan x="81.1739" y="181.925" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1075">C</tspan>
#> </text>
#> <text x="85.6259" y="175.67" id="text1076">
#>   <tspan x="85.6259" y="175.67" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1077">A</tspan>
#> </text>
#> <text x="85.0167" y="168.12" id="text1078">
#>   <tspan x="85.0167" y="168.12" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1079">C</tspan>
#> </text>
#> <text x="79.9193" y="162.353" id="text1080">
#>   <tspan x="79.9193" y="162.353" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1081">G</tspan>
#> </text>
#> <text x="79.9193" y="154.793" id="text1082">
#>   <tspan x="79.9193" y="154.793" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1083">G</tspan>
#> </text>
#> <text x="80.1293" y="147.233" id="text1084">
#>   <tspan x="80.1293" y="147.233" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1085">C</tspan>
#> </text>
#> <text x="79.9193" y="139.673" id="text1086">
#>   <tspan x="79.9193" y="139.673" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1087">G</tspan>
#> </text>
#> <text x="79.9193" y="132.113" id="text1088">
#>   <tspan x="79.9193" y="132.113" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1089">G</tspan>
#> </text>
#> <text x="86.901" y="135.474" id="text1090">
#>   <tspan x="86.901" y="135.474" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1091">U</tspan>
#> </text>
#> <text x="94.4342" y="133.611" id="text1092">
#>   <tspan x="94.4342" y="133.611" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1093">A</tspan>
#> </text>
#> <text x="98.7789" y="127.424" id="text1094">
#>   <tspan x="98.7789" y="127.424" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1095">A</tspan>
#> </text>
#> <text x="97.8374" y="119.9" id="text1096">
#>   <tspan x="97.8374" y="119.9" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1097">C</tspan>
#> </text>
#> <text x="92.5834" y="114.671" id="text1098">
#>   <tspan x="92.5834" y="114.671" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1099">A</tspan>
#> </text>
#> <text x="99.7272" y="114.671" id="text1100">
#>   <tspan x="99.7272" y="114.671" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1101">G</tspan>
#> </text>
#> <text x="107.287" y="114.671" id="text1102">
#>   <tspan x="107.287" y="114.671" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1103">G</tspan>
#> </text>
#> <text x="114.847" y="114.671" id="text1104">
#>   <tspan x="114.847" y="114.671" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1105">G</tspan>
#> </text>
#> <text x="122.407" y="114.671" id="text1106">
#>   <tspan x="122.407" y="114.671" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1107">G</tspan>
#> </text>
#> <text x="128.385" y="119.559" id="text1108">
#>   <tspan x="128.385" y="119.559" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1109">U</tspan>
#> </text>
#> <text x="135.934" y="119.962" id="text1110">
#>   <tspan x="135.934" y="119.962" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1111">U</tspan>
#> </text>
#> <text x="142.189" y="115.716" id="text1112">
#>   <tspan x="142.189" y="115.716" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1113">C</tspan>
#> </text>
#> <text x="144.392" y="108.551" id="text1114">
#>   <tspan x="144.392" y="108.551" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1115">G</tspan>
#> </text>
#> <text x="142.396" y="101.387" id="text1116">
#>   <tspan x="142.396" y="101.387" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1117">A</tspan>
#> </text>
#> <text x="136.14" y="97.1411" id="text1118">
#>   <tspan x="136.14" y="97.1411" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1119">A</tspan>
#> </text>
#> <text x="128.385" y="97.544" id="text1120">
#>   <tspan x="128.385" y="97.544" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1121">U</tspan>
#> </text>
#> <text x="122.617" y="102.431" id="text1122">
#>   <tspan x="122.617" y="102.431" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1123">C</tspan>
#> </text>
#> <text x="115.057" y="102.431" id="text1124">
#>   <tspan x="115.057" y="102.431" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1125">C</tspan>
#> </text>
#> <text x="107.497" y="102.431" id="text1126">
#>   <tspan x="107.497" y="102.431" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1127">C</tspan>
#> </text>
#> <text x="99.9372" y="102.431" id="text1128">
#>   <tspan x="99.9372" y="102.431" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1129">C</tspan>
#> </text>
#> <text x="92.3772" y="102.431" id="text1130">
#>   <tspan x="92.3772" y="102.431" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1131">U</tspan>
#> </text>
#> <text x="87.2224" y="97.101" id="text1132">
#>   <tspan x="87.2224" y="97.101" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1133">A</tspan>
#> </text>
#> <text x="86.8062" y="89.541" id="text1134">
#>   <tspan x="86.8062" y="89.541" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1135">G</tspan>
#> </text>
#> <text x="86.8062" y="81.981" id="text1136">
#>   <tspan x="86.8062" y="81.981" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1137">G</tspan>
#> </text>
#> <text x="86.8062" y="74.421" id="text1138">
#>   <tspan x="86.8062" y="74.421" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1139">G</tspan>
#> </text>
#> <text x="86.8062" y="66.861" id="text1140">
#>   <tspan x="86.8062" y="66.861" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1141">G</tspan>
#> </text>
#> <text x="87.2224" y="59.301" id="text1142">
#>   <tspan x="87.2224" y="59.301" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1143">A</tspan>
#> </text>
#> <text x="87.0162" y="51.741" id="text1144">
#>   <tspan x="87.0162" y="51.741" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1145">C</tspan>
#> </text>
#> <text x="86.8062" y="44.181" id="text1146">
#>   <tspan x="86.8062" y="44.181" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1147">G</tspan>
#> </text>
#> <text x="87.0162" y="36.621" id="text1148">
#>   <tspan x="87.0162" y="36.621" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1149">C</tspan>
#> </text>
#> <text x="87.0162" y="29.061" id="text1150">
#>   <tspan x="87.0162" y="29.061" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1151">C</tspan>
#> </text>
#> <text x="87.2224" y="21.501" id="text1152">
#>   <tspan x="87.2224" y="21.501" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1153">A</tspan>
#> </text>
#> <path fill="none" stroke="#000000" stroke-width="0.5" stroke-linecap="butt" stroke-linejoin="miter" d="M 77.4837 45.4847 A 1.78191,1.78191 270 0,0 75.7018,43.7028L 72.138 43.7028 "/>
#> <text x="65.0355" y="46.3953" id="text1154">
#>   <tspan x="65.0355" y="46.3953" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="7.1" id="tspan1155">5'</tspan>
#> </text>
#> <g id="clover-end-labels"><line x1="89.595" y1="18.935" x2="91.2384324046859" y2="4.02530081016963" stroke="black" stroke-width="0.5"/><rect x="81.4884324046859" y="-1.72469918983037" width="19.5" height="11.5" rx="3" ry="3" fill="white" stroke="black" stroke-width="0.8"/><text x="91.2384324046859" y="4.02530081016963" font-size="7.1" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="black">Glu</text></g><g id="clover-position-markers"><text x="68.2152540378444" y="111.835" font-size="5.2" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="#666666">10</text><text x="12.565" y="134.125254037844" font-size="5.2" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="#666666">20</text><text x="60.055" y="144.665" font-size="5.2" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="#666666">30</text><text x="77.295" y="168.445254037844" font-size="5.2" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="#666666">40</text><text x="107.105" y="120.765254037844" font-size="5.2" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="#666666">50</text><text x="122.094745962156" y="89.975" font-size="5.2" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="#666666">60</text><text x="99.185" y="64.295" font-size="5.2" font-family="Helvetica, Arial, sans-serif" text-anchor="middle" dominant-baseline="central" fill="#666666">70</text></g></svg></div>
# }
```
