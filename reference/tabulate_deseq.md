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


  
  {"x":{"tag":{"name":"Reactable","attribs":{"data":{"ref":["tRNA-2","tRNA-11","tRNA-3","tRNA-7","tRNA-17","tRNA-15","tRNA-4","tRNA-20","tRNA-14","tRNA-6"],"log2FoldChange":[2.64893202856241,-1.45229976886598,-1.15339838645611,0.542141549441894,-0.998432550587131,0.0402043920504232,-0.340637876031232,-0.472702482118151,-0.075574247450701,-1.27051311021181],"pvalue":[0.00108724641613662,0.00590829530265182,0.00631373389624059,0.00946227721869946,0.0199000779306516,0.0261237550061196,0.0263841792242602,0.0270977834472433,0.0275805117562413,0.0297560029197484],"padj":[0.170348634803668,0.177132067922503,0.00916986213997006,0.181237369775772,0.100249414052814,0.183330779382959,0.0762113717850298,0.177274847635999,0.129356260970235,0.184359928267077],"significant":[true,false,true,true,false,false,true,false,false,true]},"columns":[{"id":"ref","name":"ref","type":"character","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"left"},{"id":"log2FoldChange","name":"log2 FC","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["2.65","−1.45","−1.15","0.54","−1.00","0.04","−0.34","−0.47","−0.08","−1.27"],"html":true,"align":"right"},{"id":"pvalue","name":"p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["1.09&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","5.91&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","6.31&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","9.46&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.99&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.61&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.64&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.71&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.76&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.98&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"padj","name":"Adjusted p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["1.70&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.77&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","9.17&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.81&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.00&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.83&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","7.62&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.77&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.29&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.84&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>"],"html":true,"align":"right"},{"id":"significant","name":"significant","type":"logical","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"center"}],"searchable":true,"defaultPageSize":10,"showPageSizeOptions":true,"pageSizeOptions":[10,25,50,100],"paginationType":"numbers","showPagination":true,"showPageInfo":true,"minRows":1,"striped":true,"height":"auto","theme":{"color":"#333333","backgroundColor":"#FFFFFF","stripedColor":"rgba(128,128,128,0.05)","style":{"font-family":"system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif","fontSize":"16px"},"tableStyle":{"borderTopStyle":"solid","borderTopWidth":"2px","borderTopColor":"#D3D3D3"},"headerStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"groupHeaderStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"cellStyle":{"fontWeight":"normal"}},"elementId":"jobufjgfcd","dataKey":"48a62d42f8501689ccff19689e4fc105"},"children":[]},"class":"reactR_markup"},"evals":["tag.attribs.columns.0.style","tag.attribs.columns.1.style","tag.attribs.columns.2.style","tag.attribs.columns.3.style","tag.attribs.columns.4.style"],"jsHooks":[]}
```
