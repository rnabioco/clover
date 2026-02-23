# Tabulate top DESeq2 differential expression results.

Create a formatted gt table of the top significant tRNAs from
[`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
output, sorted by p-value.

## Usage

``` r
tabulate_deseq(data, lab_col = "ref", n = 10)
```

## Arguments

- data:

  A tibble from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  with at least `log2FoldChange`, `pvalue`, `padj`, and `significant`
  columns.

- lab_col:

  Column name (string) used for row labels. Default `"ref"`.

- n:

  Maximum number of rows to display. Default `10`.

## Value

An interactive `gt_tbl` object with search, sorting, and pagination.

## Examples

``` r
res <- tibble::tibble(
  ref = paste0("tRNA-", 1:20),
  log2FoldChange = rnorm(20),
  pvalue = runif(20, 0, 0.1),
  padj = runif(20, 0, 0.2),
  significant = c(rep(TRUE, 10), rep(FALSE, 10))
)
if (requireNamespace("gt", quietly = TRUE)) {
  tabulate_deseq(res)
}


  
  {"x":{"tag":{"name":"Reactable","attribs":{"data":{"ref":["tRNA-3","tRNA-5","tRNA-6","tRNA-15","tRNA-20","tRNA-10","tRNA-13","tRNA-12","tRNA-8","tRNA-14"],"log2FoldChange":[2.15875655364964,0.766983378754311,-0.308211420825515,1.4154122930517,0.480725266268846,0.322482749315297,-0.941498075931439,1.12983515341746,-0.919051596935314,0.21783764294473],"pvalue":[0.00648830933496356,0.00783842522650957,0.0133507662685588,0.0143559225834906,0.0237697737524286,0.0349695475772023,0.0360873693600297,0.0408853879896924,0.0559383239364252,0.0652510039275512],"padj":[0.0688858061097562,0.146220729034394,0.170486475294456,0.0242217760533094,0.133178157405928,0.12961655636318,0.0589250026270747,0.094627419160679,0.0226095001678914,0.122059924388304],"significant":[true,true,true,false,false,true,false,false,true,false]},"columns":[{"id":"ref","name":"ref","type":"character","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"left"},{"id":"log2FoldChange","name":"log2 FC","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["2.16","0.77","−0.31","1.42","0.48","0.32","−0.94","1.13","−0.92","0.22"],"html":true,"align":"right"},{"id":"pvalue","name":"p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["6.49&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","7.84&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.34&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.44&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.38&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.50&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.61&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.09&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.59&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","6.53&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"padj","name":"Adjusted p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["6.89&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.46&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.70&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","2.42&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.33&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.30&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","5.89&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","9.46&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.26&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.22&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>"],"html":true,"align":"right"},{"id":"significant","name":"significant","type":"logical","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"center"}],"searchable":true,"defaultPageSize":10,"showPageSizeOptions":true,"pageSizeOptions":[10,25,50,100],"paginationType":"numbers","showPagination":true,"showPageInfo":true,"minRows":1,"striped":true,"height":"auto","theme":{"color":"#333333","backgroundColor":"#FFFFFF","stripedColor":"rgba(128,128,128,0.05)","style":{"font-family":"system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif","fontSize":"16px"},"tableStyle":{"borderTopStyle":"solid","borderTopWidth":"2px","borderTopColor":"#D3D3D3"},"headerStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"groupHeaderStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"cellStyle":{"fontWeight":"normal"}},"elementId":"imnwerlwzl","dataKey":"6a896986ab2e4e1aa0a623450489bd41"},"children":[]},"class":"reactR_markup"},"evals":["tag.attribs.columns.0.style","tag.attribs.columns.1.style","tag.attribs.columns.2.style","tag.attribs.columns.3.style","tag.attribs.columns.4.style"],"jsHooks":[]}
```
