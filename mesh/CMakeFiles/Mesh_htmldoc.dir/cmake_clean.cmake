file(REMOVE_RECURSE
  "mesh_kernel.inc"
  "reduce_kernel.inc"
  "install_manifest.txt"
  "docs/_build"
  "docs/htmldoc.out"
  "docs/pdfdoc.out"
  "docs/singlehtmldoc.out"
  "CMakeFiles/Mesh_htmldoc"
  "_build/html/UsersGuide.html"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/Mesh_htmldoc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
