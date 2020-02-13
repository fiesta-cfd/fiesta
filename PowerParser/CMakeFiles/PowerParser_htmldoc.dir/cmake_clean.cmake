file(REMOVE_RECURSE
  "FParser_module.mod"
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/PowerParser_htmldoc"
  "_build/html/UsersGuide.html"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/PowerParser_htmldoc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
