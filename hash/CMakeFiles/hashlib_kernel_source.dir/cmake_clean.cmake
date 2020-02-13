file(REMOVE_RECURSE
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/hashlib_kernel_source"
  "hashlib_kern.inc"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/hashlib_kernel_source.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
