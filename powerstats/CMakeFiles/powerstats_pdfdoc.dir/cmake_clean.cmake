file(REMOVE_RECURSE
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/powerstats_pdfdoc"
  "_build/latex/powerstats.pdf"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/powerstats_pdfdoc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
