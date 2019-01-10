if(TARGET colors::colors)
	return()
endif()

set(MY_FLAGS
		-fdiagnostics-color=auto
)

# Flags above don't make sense for MSVC
if(MSVC)
	set(MY_FLAGS)
endif()

include(CheckCXXCompilerFlag)

add_library(colors_colors INTERFACE)
add_library(colors::colors ALIAS colors_colors)

foreach(FLAG IN ITEMS ${MY_FLAGS})
	string(REPLACE "=" "-" FLAG_VAR "${FLAG}")
	if(NOT DEFINED IS_SUPPORTED_${FLAG_VAR})
		check_cxx_compiler_flag("${FLAG}" IS_SUPPORTED_${FLAG_VAR})
	endif()
	if(IS_SUPPORTED_${FLAG_VAR})
		target_compile_options(colors_colors INTERFACE ${FLAG})
	endif()
endforeach()
