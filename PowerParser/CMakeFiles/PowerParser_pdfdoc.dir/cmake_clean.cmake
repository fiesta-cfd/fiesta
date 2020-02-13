file(REMOVE_RECURSE
  "FParser_module.mod"
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/PowerParser_pdfdoc"
  "_build/latex/PowerParser.pdf"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/PowerParser_pdfdoc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
