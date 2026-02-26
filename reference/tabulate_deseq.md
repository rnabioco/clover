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


  
  {"x":{"tag":{"name":"Reactable","attribs":{"data":{"ref":["tRNA-10","tRNA-2","tRNA-6","tRNA-20","tRNA-4","tRNA-18","tRNA-19","tRNA-13","tRNA-3","tRNA-14"],"log2FoldChange":[0.149405696123499,-2.37992516827699,0.0361095616099068,-1.21047936165691,0.127840486753431,0.709065963534284,1.04692267028002,0.392105946790304,0.597565620228058,0.469933174095712],"pvalue":[0.00344292330555618,0.0056073043262586,0.0121108880266547,0.0145407346775755,0.0294625013135374,0.0398619815241545,0.043052065372467,0.0451858439715579,0.0473137095803395,0.0516804101876915],"padj":[0.105741854477674,0.115756298322231,0.126908854441717,0.128136096987873,0.148465236043558,0.18366756895557,0.0297420693561435,0.0218097571283579,0.0140834241174161,0.0139891034457833],"significant":[true,true,true,false,true,false,false,false,true,false]},"columns":[{"id":"ref","name":"ref","type":"character","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"left"},{"id":"log2FoldChange","name":"log2 FC","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["0.15","−2.38","0.04","−1.21","0.13","0.71","1.05","0.39","0.60","0.47"],"html":true,"align":"right"},{"id":"pvalue","name":"p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["3.44&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","5.61&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.21&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.45&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.95&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.99&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.31&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.52&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.73&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.17&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"padj","name":"Adjusted p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["1.06&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.16&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.27&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.28&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.48&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.84&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","2.97&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.18&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.41&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.40&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"significant","name":"significant","type":"logical","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"center"}],"searchable":true,"defaultPageSize":10,"showPageSizeOptions":true,"pageSizeOptions":[10,25,50,100],"paginationType":"numbers","showPagination":true,"showPageInfo":true,"minRows":1,"striped":true,"height":"auto","theme":{"color":"#333333","backgroundColor":"#FFFFFF","stripedColor":"rgba(128,128,128,0.05)","style":{"font-family":"system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif","fontSize":"16px"},"tableStyle":{"borderTopStyle":"solid","borderTopWidth":"2px","borderTopColor":"#D3D3D3"},"headerStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"groupHeaderStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"cellStyle":{"fontWeight":"normal"}},"elementId":"ehmulxtlvm","dataKey":"214a4479abe881af6d72a66319aeb7f5"},"children":[]},"class":"reactR_markup"},"evals":["tag.attribs.columns.0.style","tag.attribs.columns.1.style","tag.attribs.columns.2.style","tag.attribs.columns.3.style","tag.attribs.columns.4.style"],"jsHooks":[]}
```
