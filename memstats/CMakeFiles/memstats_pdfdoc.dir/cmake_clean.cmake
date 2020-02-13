file(REMOVE_RECURSE
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/memstats_pdfdoc"
  "_build/latex/memstats.pdf"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/memstats_pdfdoc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
