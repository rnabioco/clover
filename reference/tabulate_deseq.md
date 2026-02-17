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


  
  {"x":{"tag":{"name":"Reactable","attribs":{"data":{"ref":["tRNA-9","tRNA-2","tRNA-16","tRNA-1","tRNA-14","tRNA-11","tRNA-3","tRNA-8","tRNA-15","tRNA-7"],"log2FoldChange":[-0.8286035327893,0.441428225928461,-0.010303318628592,-0.0778868243875754,0.1496793540957,-0.48005348427113,0.128922896223691,1.45584140300823,-1.4333211002982,-0.751723322890923],"pvalue":[0.00147188026458025,0.00581756506580859,0.00746336323209107,0.0166490767849609,0.0360042401123792,0.0453089235117659,0.0459139019018039,0.0496825250331312,0.0503146104514599,0.0524280134355649],"padj":[0.0358068581204861,0.112877145409584,0.14169685235247,0.152211881009862,0.165584264416248,0.142682390334085,0.155690860701725,0.0913498361594975,0.128613806329668,0.168846283247694],"significant":[true,true,false,true,false,false,true,true,false,true]},"columns":[{"id":"ref","name":"ref","type":"character","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"left"},{"id":"log2FoldChange","name":"log2 FC","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["−0.83","0.44","−0.01","−0.08","0.15","−0.48","0.13","1.46","−1.43","−0.75"],"html":true,"align":"right"},{"id":"pvalue","name":"p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["1.47&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","5.82&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","7.46&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−3<\/sup>","1.66&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","3.60&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.53&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.59&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","4.97&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.03&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","5.24&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>"],"html":true,"align":"right"},{"id":"padj","name":"Adjusted p-value","type":"numeric","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","cell":["3.58&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.13&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.42&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.52&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.66&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.43&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.56&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","9.13&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−2<\/sup>","1.29&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>","1.69&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−1<\/sup>"],"html":true,"align":"right"},{"id":"significant","name":"significant","type":"logical","na":"NA","minWidth":125,"style":"function(rowInfo, colInfo) {\nconst rowIndex = rowInfo.index + 1\n}","html":true,"align":"center"}],"searchable":true,"defaultPageSize":10,"showPageSizeOptions":true,"pageSizeOptions":[10,25,50,100],"paginationType":"numbers","showPagination":true,"showPageInfo":true,"minRows":1,"striped":true,"height":"auto","theme":{"color":"#333333","backgroundColor":"#FFFFFF","stripedColor":"rgba(128,128,128,0.05)","style":{"font-family":"system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif","fontSize":"16px"},"tableStyle":{"borderTopStyle":"solid","borderTopWidth":"2px","borderTopColor":"#D3D3D3"},"headerStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"groupHeaderStyle":{"fontWeight":"normal","backgroundColor":"transparent","borderBottomStyle":"solid","borderBottomWidth":"2px","borderBottomColor":"#D3D3D3"},"cellStyle":{"fontWeight":"normal"}},"elementId":"sijyaxpvjf","dataKey":"9c1f49faf9a23812347f10be38edf412"},"children":[]},"class":"reactR_markup"},"evals":["tag.attribs.columns.0.style","tag.attribs.columns.1.style","tag.attribs.columns.2.style","tag.attribs.columns.3.style","tag.attribs.columns.4.style"],"jsHooks":[]}
```
