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


  
  {"x":{"tag":{"name":"Reactable","attribs":{"data":{"ref":["tRNA-7","tRNA-20","tRNA-17","tRNA-5","tRNA-16","tRNA-2","tRNA-19","tRNA-6","tRNA-11","tRNA-4"],"log2FoldChange":[-1.2741013109596,-0.00795800684151624,-1.57027435202884,-0.652674974433022,0.0883923627751587,0.569298098772819,0.652678559530207,-0.147169331518372,-0.820979907165238,1.24875501630229],"pvalue":[0.00746336323209107,0.0179034290602431,0.0359530670568347,0.0360042401123792,0.0378960735630244,0.0453089235117659,0.0456749180797488,0.0503146104514599,0.0536133256973699,0.0544705166015774],"padj":[0.14169685235247,0.0824512774124742,0.00173161495476961,0.165584264416248,0.0701176324859262,0.142682390334085,0.144987019430846,0.128613806329668,0.0219340038020164,0.125291298842058],"significant":[true,false,false,true,false,true,false,true,false,true]},"columns":[{"id":"ref","name":"ref","type":"character","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"left"},{"id":"log2FoldChange","name":"log2 FC","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["−1.27","−0.01","−1.57","−0.65","0.09","0.57","0.65","−0.15","−0.82","1.25"],"html":true,"align":"right"},{"id":"pvalue","name":"p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["7.46&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.79&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.60&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.60&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.79&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.53&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.57&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.03&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.36&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.45&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"padj","name":"Adjusted p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["1.42&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","8.25&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.73&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.66&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","7.01&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.43&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.45&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.29&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","2.19&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.25&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>"],"html":true,"align":"right"},{"id":"significant","name":"significant","type":"logical","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"center"}],"searchable":true,"defaultPageSize":10,"showPageSizeOptions":true,"pageSizeOptions":[10,25,50,100],"paginationType":"numbers","showPagination":true,"showPageInfo":true,"minRows":1,"striped":true,"height":"auto","theme":{"color":"#333333","backgroundColor":"#FFFFFF","stripedColor":"rgba(128,128,128,0.05)","style":{"font-family":"system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif","fontSize":"16px"},"tableStyle":{"borderTopStyle":"solid","borderTopWidth":"2px","borderTopColor":"#D3D3D3"},"headerStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"groupHeaderStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"cellStyle":{"fontWeight":"normal"}},"elementId":"fhcepyvzqn","dataKey":"2a6dbe83a47d588cec68cfdcb357a270"},"children":[]},"class":"reactR_markup"},"evals":["tag.attribs.columns.0.style","tag.attribs.columns.1.style","tag.attribs.columns.2.style","tag.attribs.columns.3.style","tag.attribs.columns.4.style"],"jsHooks":[]}
```
