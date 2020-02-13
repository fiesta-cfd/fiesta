file(REMOVE_RECURSE
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/Crux_doc"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/Crux_doc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
