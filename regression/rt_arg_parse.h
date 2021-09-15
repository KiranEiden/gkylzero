#pragma once

#include <gkyl_basis.h>

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define APP_ARGS_DEFAULT_FILE_NAME "name-not-set"

struct gkyl_app_args {
  bool use_gpu; // should this be run on GPU?
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  int num_steps; // number of steps
  char file_name[1024]; // name of input file
  enum gkyl_basis_type basis_type; // type of basis functions to use
};

static int
get_basis_type(const char *nm)
{
  if (strcmp(nm, "ms") == 0) {
    return GKYL_BASIS_MODAL_SERENDIPITY;
  }
  else if (strcmp(nm, "mt") == 0) {
    return GKYL_BASIS_MODAL_TENSOR;
  }
  return -1;
}

static struct gkyl_app_args
parse_app_args(int argc, char **argv)
{
  bool use_gpu = false;
  bool step_mode = false;
  int num_steps = INT_MAX;

  struct gkyl_app_args args;

  strcpy(args.file_name, APP_ARGS_DEFAULT_FILE_NAME); // default
  args.basis_type = GKYL_BASIS_MODAL_SERENDIPITY;

  int c;
  while ((c = getopt(argc, argv, "+hgs:i:b:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: <app_name> -g -s num_steps -i inp_file -b [ms|mt]\n");
        printf(" -g     Run on GPUs if GPUs are present and code built for GPUs\n");
        printf(" -sN    Only run N steps of simulation\n");
        printf(" -b     Basis function to use (ms: Modal serendipity; mt: Modal tensor-product)\n");
        printf("        (Ignored for finite-volume solvers)\n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
        break;        
      
      case 's':
        step_mode = true;
        num_steps = atoi(optarg);
        break;

      case 'i':
        strcpy(args.file_name, optarg);
        break;

      case 'b':
        args.basis_type = get_basis_type(optarg);
        assert(args.basis_type != -1);
        break;

      case '?':
        break;
    }
  }
  
  args.use_gpu = use_gpu;
  args.step_mode = step_mode;
  args.num_steps = num_steps;

  return args;
}
