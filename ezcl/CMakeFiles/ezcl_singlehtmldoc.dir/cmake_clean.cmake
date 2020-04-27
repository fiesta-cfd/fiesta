file(REMOVE_RECURSE
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/ezcl_singlehtmldoc"
  "_build/singlehtml/index.html"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/ezcl_singlehtmldoc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
