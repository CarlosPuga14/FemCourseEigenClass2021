
add_library(course_library  ${SOURCES} ${HEADERS})
target_link_libraries(course_library PRIVATE Eigen3::Eigen)
# Mac OS X specific code
if (APPLE)
  target_compile_definitions(course_library PRIVATE MACOSX)
endif()


add_executable(TestEigen mainEigenTest.cpp)
target_link_libraries (TestEigen course_library Eigen3::Eigen)

add_executable(Fem2021 main2021.cpp)
target_link_libraries (Fem2021 Eigen3::Eigen course_library)

add_executable(Contribute mainContribute.cpp)
target_link_libraries (Contribute Eigen3::Eigen course_library)

add_executable(PostProcess mainPostProcess.cpp)
target_link_libraries (PostProcess Eigen3::Eigen course_library)

add_executable(Calcstiff mainCalcstiff.cpp)
target_link_libraries (Calcstiff Eigen3::Eigen course_library)

add_executable(Assemble mainAssemble.cpp)
target_link_libraries (Assemble Eigen3::Eigen course_library)

add_executable(TestIntegrate TestIntegrate.cpp)
target_link_libraries (TestIntegrate  PRIVATE course_library Eigen3::Eigen Catch2::Catch2)

add_executable(TestGeom TestGeom.cpp)
target_link_libraries (TestGeom  PRIVATE course_library Eigen3::Eigen Catch2::Catch2)

add_executable(TestShape TestShape.cpp)
target_link_libraries (TestShape PRIVATE course_library Eigen3::Eigen Catch2::Catch2)

add_executable(TestOneDProblem TestOneDProblem.cpp)
target_link_libraries (TestOneDProblem  course_library Eigen3::Eigen)

add_executable(mainCont mainContribute.cpp)
target_link_libraries (mainCont course_library Eigen3::Eigen)

add_executable(mainPostPos mainPosProcessamento.cpp)
target_link_libraries (mainPostPos course_library Eigen3::Eigen)

add_executable(TestTwoDProblem TestTwoDProblem.cpp)
target_link_libraries (TestTwoDProblem course_library Eigen3::Eigen)

add_executable(TestThreeDProblem TestThreeDProblem.cpp)
target_link_libraries (TestThreeDProblem course_library Eigen3::Eigen)

# Mac OS X specific code
if (APPLE)
  target_compile_definitions(TestTwoDProblem PRIVATE MACOSX)
  target_compile_definitions(TestOneDProblem PRIVATE MACOSX)
  target_compile_definitions(TestThreeDProblem PRIVATE MACOSX)
endif()



file(GLOB MY_MESH_FILES
  "*.msh" "*.geo"
)

file(COPY ${MY_MESH_FILES} DESTINATION .)

# configure_file(*.msh *.msh COPYONLY)
#
# configure_file(*.geo *.geo COPYONLY)
