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


  
  {"x":{"tag":{"name":"Reactable","attribs":{"data":{"ref":["tRNA-14","tRNA-11","tRNA-20","tRNA-8","tRNA-7","tRNA-10","tRNA-2","tRNA-4","tRNA-15","tRNA-6"],"log2FoldChange":[0.239959572121817,1.893360463671,1.87774387154249,-0.21223603483794,-0.010303318628592,-2.10215247889527,-0.48005348427113,1.46011018039623,0.0608988932499574,-1.4333211002982],"pvalue":[0.00248860311694443,0.0179034290602431,0.0281812581699342,0.0359530670568347,0.0378960735630244,0.0456749180797488,0.0536133256973699,0.056438572704792,0.0626456494210288,0.0639060988556594],"padj":[0.00860025864094496,0.0824512774124742,0.15242218291387,0.00173161495476961,0.0701176324859262,0.144987019430846,0.0219340038020164,0.170144265471026,0.13692882261239,0.135132380388677],"significant":[false,false,false,true,true,true,true,true,false,true]},"columns":[{"id":"ref","name":"ref","type":"character","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"left"},{"id":"log2FoldChange","name":"log2 FC","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["0.24","1.89","1.88","−0.21","−0.01","−2.10","−0.48","1.46","0.06","−1.43"],"html":true,"align":"right"},{"id":"pvalue","name":"p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["2.49&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.79&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","2.82&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.60&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.79&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.57&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.36&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.64&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","6.26&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","6.39&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"padj","name":"Adjusted p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["8.60&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","8.25&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.52&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.73&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","7.01&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.45&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","2.19&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.70&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.37&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.35&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>"],"html":true,"align":"right"},{"id":"significant","name":"significant","type":"logical","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"center"}],"searchable":true,"defaultPageSize":10,"showPageSizeOptions":true,"pageSizeOptions":[10,25,50,100],"paginationType":"numbers","showPagination":true,"showPageInfo":true,"minRows":1,"striped":true,"height":"auto","theme":{"color":"#333333","backgroundColor":"#FFFFFF","stripedColor":"rgba(128,128,128,0.05)","style":{"font-family":"system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif","fontSize":"16px"},"tableStyle":{"borderTopStyle":"solid","borderTopWidth":"2px","borderTopColor":"#D3D3D3"},"headerStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"groupHeaderStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"cellStyle":{"fontWeight":"normal"}},"elementId":"qnitvfqkck","dataKey":"f10de0fdccc8da755c91a6e44b544fae"},"children":[]},"class":"reactR_markup"},"evals":["tag.attribs.columns.0.style","tag.attribs.columns.1.style","tag.attribs.columns.2.style","tag.attribs.columns.3.style","tag.attribs.columns.4.style"],"jsHooks":[]}
```
