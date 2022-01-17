## Standing comments
- We occasionally see a NOTE about "(possibly) invalid URLs," but all of the URLs/DOIs we reference are valid.


## CRAN build NOTE from January 17, 2022
Version: 2.6.0
Check: dependencies in R code
Result: NOTE
    Namespace in Imports field not imported from: ‘mathjaxr’
     All declared Imports should be used.
Flavor: r-oldrel-macos-x86_64

Comment: `mathjaxr` is needed to properly render our documentation.
