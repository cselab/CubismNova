
# parser
test_data = custom_target('ini-parser-data',
  command: ['rsync', '@INPUT@', '@OUTPUT@'],
  input: files([
    'INIParser/INIFiles/incA.ini',
    'INIParser/INIFiles/incAa.ini',
    'INIParser/INIFiles/incAaa.ini',
    'INIParser/INIFiles/incB.ini',
    'INIParser/INIFiles/main.ini',
    ]),
  output: 'INIFiles',
)
p = executable('ini-parser',
  [files([
    'INIParser/INIParserTest.cpp',
    ]), test_data],
  include_directories: inc_cubism_nova,
  dependencies: gtest_main_dep,
  link_with: lib_cubism_util,
)
test('ini-parser', p,
  workdir: meson.current_build_dir(),
  protocol: 'gtest',
  suite: 'Util',
)

# version
v = executable('version',
  files([
    'Version/VersionTest.cpp',
    ]),
  include_directories: inc_cubism_nova,
  dependencies: gtest_main_dep,
  link_with: lib_cubism_util,
)
test('version', v,
  protocol: 'gtest',
  suite: 'Util',
)
