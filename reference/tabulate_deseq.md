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


  
  {"x":{"tag":{"name":"Reactable","attribs":{"data":{"ref":["tRNA-7","tRNA-10","tRNA-13","tRNA-6","tRNA-16","tRNA-12","tRNA-5","tRNA-17","tRNA-8","tRNA-20"],"log2FoldChange":[0.675270015431045,0.372003230867959,-0.302706284526524,-0.248198556052917,-2.21814114794348,1.0439388958914,0.710346515685883,-0.275302017145566,-0.595348959958748,0.258902502532157],"pvalue":[0.00683417033869773,0.0116309177130461,0.0133348547620699,0.0203208829043433,0.0203664392232895,0.0226108731469139,0.024249267578125,0.0256982919061556,0.0307273737154901,0.0315594685263932],"padj":[0.17083473042585,0.0991779293399304,0.0878299999050796,0.11189823942259,0.00355399120599031,0.0831922901328653,0.00532493400387466,0.145635994663462,0.015176609903574,0.0332981535699218],"significant":[true,true,false,true,false,false,true,false,true,false]},"columns":[{"id":"ref","name":"ref","type":"character","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"left"},{"id":"log2FoldChange","name":"log2 FC","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["0.68","0.37","−0.30","−0.25","−2.22","1.04","0.71","−0.28","−0.60","0.26"],"html":true,"align":"right"},{"id":"pvalue","name":"p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["6.83&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.16&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.33&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.03&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.04&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.26&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.42&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.57&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.07&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.16&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"padj","name":"Adjusted p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["1.71&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","9.92&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","8.78&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.12&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","3.55&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","8.32&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.32&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.46&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.52&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.33&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"significant","name":"significant","type":"logical","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"center"}],"searchable":true,"defaultPageSize":10,"showPageSizeOptions":true,"pageSizeOptions":[10,25,50,100],"paginationType":"numbers","showPagination":true,"showPageInfo":true,"minRows":1,"striped":true,"height":"auto","theme":{"color":"#333333","backgroundColor":"#FFFFFF","stripedColor":"rgba(128,128,128,0.05)","style":{"font-family":"system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif","fontSize":"16px"},"tableStyle":{"borderTopStyle":"solid","borderTopWidth":"2px","borderTopColor":"#D3D3D3"},"headerStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"groupHeaderStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"cellStyle":{"fontWeight":"normal"}},"elementId":"ekpgwxpefr","dataKey":"73f9e9b21b7fecf0d9a1d93abdf46c71"},"children":[]},"class":"reactR_markup"},"evals":["tag.attribs.columns.0.style","tag.attribs.columns.1.style","tag.attribs.columns.2.style","tag.attribs.columns.3.style","tag.attribs.columns.4.style"],"jsHooks":[]}
```
