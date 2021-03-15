# vcs_dep = vcs_tag(
#   command: ['python3', join_paths(meson.project_source_root(),
#   'tools', 'version', 'gen_version.py'), '--project_root',
#   meson.project_source_root(), '--input', '@INPUT@', '--output',
#   '@OUTPUT@'],
#   input: 'Version.cpp.in',
#   output: 'Version.cpp'
# )
vcs_dep = vcs_tag(
  input: 'Version.cpp.in',
  output: 'Version.cpp'
)
