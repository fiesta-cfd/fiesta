file(REMOVE_RECURSE
  "l7_kernel.inc"
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/L7_htmldoc"
  "_build/html/UsersGuide.html"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/L7_htmldoc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
