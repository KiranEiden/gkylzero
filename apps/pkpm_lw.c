#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_lua_utils.h>
#include <gkyl_pkpm.h>
#include <gkyl_pkpm_lw.h>
#include <gkyl_pkpm_priv.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

// Check and fetch user-data based on metatable name
#define CHECK_UDATA(L, mnm) luaL_checkudata(L, 1, mnm)

// For debugging
#define trace_stack_top(L, fnm) do { \
      fprintf(stdout, "Inside function %s\n", fnm);                              \
      fprintf(stdout, "--> Top of stack is %s\n", lua_typename(L, lua_type(L, -1))); \
    } while (0);

// Get basis type from string
static enum gkyl_basis_type
get_basis_type(const char *bnm)
{
  if (strcmp(bnm, "serendipity") == 0)
    return GKYL_BASIS_MODAL_SERENDIPITY;
  if (strcmp(bnm, "tensor") == 0)
    return GKYL_BASIS_MODAL_TENSOR;
  if (strcmp(bnm, "hybrid") == 0)
    return GKYL_BASIS_MODAL_HYBRID;

  return GKYL_BASIS_MODAL_SERENDIPITY;
}

// Magic IDs for use in distinguishing various species and field types
enum pkpm_magic_ids {
  PKPM_SPECIES_DEFAULT = 100, // non-relativistic PKPM model
  PKPM_FIELD_DEFAULT, // Maxwell equations
  PKPM_COLLISIONS_DEFAULT, // LBO Collisions
};

// Used in call back passed to the initial conditions
struct lua_func_ctx {
  int func_ref; // reference to Lua function in registery
  int ndim, nret; // dimensions of function, number of return values
  lua_State *L; // Lua state
};

static void
eval_ic(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct lua_func_ctx *fr = ctx;
  lua_State *L = fr->L;

  int ndim = fr->ndim;
  int nret = fr->nret;
  lua_rawgeti(L, LUA_REGISTRYINDEX, fr->func_ref);
  lua_pushnumber(L, t);
  lua_createtable(L, GKYL_MAX_DIM, 0);

  for (int i=0; i<ndim; ++i) {
    lua_pushnumber(L, xn[i]);
    lua_rawseti(L, -2, i+1); 
  }

  if (lua_pcall(L, 2, nret, 0)) {
    const char* ret = lua_tostring(L, -1);
    luaL_error(L, "*** eval_ic ERROR: %s\n", ret);
  }

  for (int i=nret-1; i>=0; --i) { // need to fetch in reverse order
    fout[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
}

/* *****************/
/* Collisions methods */
/* *****************/

// Metatable name for field input struct
#define PKPM_COLLISIONS_METATABLE_NM "GkeyllZero.PKPM.Collisions"

// Lua userdata object for constructing field input
struct pkpm_collisions_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_pkpm_collisions pkpm_collisions; // input struct to construct collisions
  struct lua_func_ctx init_ref; // Lua registery reference to initilization function for (self) collision frequency
};

static int
pkpm_collisions_lw_new(lua_State *L)
{
  printf("Inside pkpm_collisions_lw_new\n");
  struct gkyl_pkpm_collisions pkpm_collisions = { };

  pkpm_collisions.collision_id = GKYL_LBO_COLLISIONS;  

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "self_nu"))
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Collisions must have an \"self_nu\" function for collision frequency!");

  struct pkpm_collisions_lw *pkpmc_lw = lua_newuserdata(L, sizeof(*pkpmc_lw));

  pkpmc_lw->magic = PKPM_COLLISIONS_DEFAULT;
  pkpmc_lw->pkpm_collisions = pkpm_collisions;
  
  pkpmc_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };  
  
  // set metatable
  luaL_getmetatable(L, PKPM_COLLISIONS_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg pkpm_collisions_ctor[] = {
  { "new",  pkpm_collisions_lw_new },
  { 0, 0 }
};

/* *****************/
/* Species methods */
/* *****************/

// Metatable name for species input struct
#define PKPM_SPECIES_METATABLE_NM "GkeyllZero.PKPM.Species"

// Lua userdata object for constructing species input
struct pkpm_species_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_pkpm_species pkpm_species; // input struct to construct species
  int vdim; // velocity dimensions
  bool evolve; // is this species evolved?
  struct lua_func_ctx init_dist_ref; // Lua registery reference to initilization distribution function
  struct lua_func_ctx init_fluid_ref; // Lua registery reference to initilization momentum

  struct pkpm_collisions_lw *collisions_lw; // pointer to Lua collisions table
};

static int
pkpm_species_lw_new(lua_State *L)
{
  printf("Inside pkpm_species_lw_new\n");
  int vdim  = 0;
  struct gkyl_pkpm_species pkpm_species = { };

  pkpm_species.model_id = GKYL_MODEL_DEFAULT;
  
  pkpm_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  pkpm_species.mass = glua_tbl_get_number(L, "mass", 1.0);

  with_lua_tbl_tbl(L, "cells") {
    vdim = glua_objlen(L);
    for (int d=0; d<vdim; ++d)
      pkpm_species.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d=0; d<vdim; ++d)
      pkpm_species.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<vdim; ++d)
      pkpm_species.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  int init_dist_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init_dist"))
    init_dist_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Species must have an \"init_dist\" function for initial conditions!");

  int init_fluid_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init_fluid"))
    init_fluid_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Species must have an \"init_fluid\" function for initial conditions!");

  struct pkpm_collisions_lw *pkpmc = 0;
  with_lua_tbl_key(L, "collisions") {
    printf(".. Reading collisions ....\n");
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct pkpm_collisions_lw *my_pkpmc = lua_touserdata(L, -1);
      if (my_pkpmc->magic == PKPM_COLLISIONS_DEFAULT)
        pkpmc = my_pkpmc;
    }
  }  
  
  struct pkpm_species_lw *pkpms_lw = lua_newuserdata(L, sizeof(*pkpms_lw));
  pkpms_lw->magic = PKPM_SPECIES_DEFAULT;
  pkpms_lw->vdim = vdim;
  pkpms_lw->evolve = evolve;
  pkpms_lw->pkpm_species = pkpm_species;

  pkpms_lw->init_dist_ref = (struct lua_func_ctx) {
    .func_ref = init_dist_ref,
    .ndim = 0, // this will be set later
    .nret = 2,
    .L = L,
  };

  pkpms_lw->init_fluid_ref = (struct lua_func_ctx) {
    .func_ref = init_fluid_ref,
    .ndim = 0, // this will be set later
    .nret = 3,
    .L = L,
  };

  pkpms_lw->collisions_lw = pkpmc;
   
  // set metatable
  luaL_getmetatable(L, PKPM_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg pkpm_species_ctor[] = {
  {"new", pkpm_species_lw_new},
  {0, 0}
};

/* *****************/
/* Field methods */
/* *****************/

// Metatable name for field input struct
#define PKPM_FIELD_METATABLE_NM "GkeyllZero.PKPM.Field"

// Lua userdata object for constructing field input
struct pkpm_field_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_pkpm_field pkpm_field; // input struct to construct field
  bool evolve; // is this field evolved?
  struct lua_func_ctx init_ref; // Lua registery reference to initilization function
};

static int
pkpm_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_pkpm_field pkpm_field = { };

  pkpm_field.field_id = GKYL_FIELD_E_B;  
  
  pkpm_field.epsilon0 = glua_tbl_get_number(L, "epsilon0", 1.0);
  pkpm_field.mu0 = glua_tbl_get_number(L, "mu0", 1.0);
  pkpm_field.elcErrorSpeedFactor = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  pkpm_field.mgnErrorSpeedFactor = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init"))
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Field must have an \"init\" function for initial conditions!");

  struct pkpm_field_lw *pkpmf_lw = lua_newuserdata(L, sizeof(*pkpmf_lw));

  pkpmf_lw->magic = PKPM_FIELD_DEFAULT;
  pkpmf_lw->evolve = evolve;
  pkpmf_lw->pkpm_field = pkpm_field;
  
  pkpmf_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // this will be set later
    .nret = 6,
    .L = L,
  };  
  
  // set metatable
  luaL_getmetatable(L, PKPM_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg pkpm_field_ctor[] = {
  { "new",  pkpm_field_lw_new },
  { 0, 0 }
};

/* *************/
/* App methods */
/* *************/

// Metatable name for top-level PKPM App
#define PKPM_APP_METATABLE_NM "GkeyllZero.PKPM.App"

// Lua userdata object for holding PKPM app and run parameters
struct pkpm_app_lw {
  gkyl_pkpm_app *app; // PKPM app object
  struct lua_func_ctx species_dist_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' distribution function
  struct lua_func_ctx species_fluid_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' momentum
  struct lua_func_ctx collisions_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' collision frequency
  struct lua_func_ctx field_func_ctx; // function context for field
  
  double tstart, tend; // start and end times of simulation
  int nframe; // number of data frames to write
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct pkpm_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1};
  
  int curr = 0;
  lua_pushnil(L); // initial key is nil
  while (lua_next(L, TKEY) != 0) {
    // key at TKEY and value at TVAL
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct pkpm_species_lw *pkpms = lua_touserdata(L, TVAL);
      if (pkpms->magic == PKPM_SPECIES_DEFAULT) {
        
        pkpms->init_dist_ref.ndim = cdim + pkpms->vdim;
        pkpms->init_fluid_ref.ndim = cdim;
        pkpms->collisions_lw->init_ref.ndim = cdim;
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(pkpms->pkpm_species.name, key);
        }
        species[curr++] = pkpms;
      }
    }
    lua_pop(L, 1);
  }
  return curr;
}

// Create top-level App object
static int
pkpm_app_new(lua_State *L)
{
  struct pkpm_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-pkpm" is used.
  const char *sim_name = "g0-pkpm";
  with_lua_global(L, "GKYL_OUT_PREFIX") {
    if (lua_isstring(L, -1))
      sim_name = lua_tostring(L, -1);
  }
  
  // initialize app using table inputs (table is on top of stack)

  app_lw->tstart = glua_tbl_get_number(L, "tStart", 0.0);
  app_lw->tend = glua_tbl_get_number(L, "tEnd", 1.0);
  app_lw->nframe = glua_tbl_get_integer(L, "nFrame", 1);

  struct gkyl_pkpm pkpm = { }; // input table for app

  strcpy(pkpm.name, sim_name);
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    pkpm.cdim = cdim = glua_objlen(L);
    for (int d=0; d<cdim; ++d)
      pkpm.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d=0; d<cdim; ++d)
      pkpm.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<cdim; ++d)
      pkpm.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  pkpm.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  pkpm.basis_type = get_basis_type(
    glua_tbl_get_string(L, "basis", "serendipity")
  );

  pkpm.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      pkpm.num_periodic_dir = glua_objlen(L);
      for (int d=0; d<pkpm.num_periodic_dir; ++d)
        // indexes are off by 1 between Lua and C
        pkpm.periodic_dirs[d] = glua_tbl_iget_integer(L, d+1, 0)-1;
    }
  }

  struct pkpm_species_lw *species[GKYL_MAX_SPECIES];
  // set all species input
  pkpm.num_species = get_species_inp(L, cdim, species);
  for (int s=0; s<pkpm.num_species; ++s) {
    pkpm.species[s] = species[s]->pkpm_species;
    pkpm.vdim = species[s]->vdim;
    
    // get context from distribution functions and momentum
    app_lw->species_dist_func_ctx[s] = species[s]->init_dist_ref;
    app_lw->species_fluid_func_ctx[s] = species[s]->init_fluid_ref;
    pkpm.species[s].init_dist = eval_ic;
    pkpm.species[s].init_fluid = eval_ic;
    pkpm.species[s].ctx_dist = &app_lw->species_dist_func_ctx[s];
    pkpm.species[s].ctx_fluid = &app_lw->species_fluid_func_ctx[s];

    pkpm.species[s].collisions = species[s]->collisions_lw->pkpm_collisions;

    // get context for (self) collision frequency
    app_lw->collisions_func_ctx[s] = species[s]->collisions_lw->init_ref;
    pkpm.species[s].collisions.self_nu = eval_ic;
    pkpm.species[s].collisions.ctx = &app_lw->collisions_func_ctx[s];
  }

  // set field input
  pkpm.skip_field = true;
  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct pkpm_field_lw *pkpmf = lua_touserdata(L, -1);
      if (pkpmf->magic == PKPM_FIELD_DEFAULT) {
        pkpmf->init_ref.ndim = cdim;

        pkpm.field = pkpmf->pkpm_field;
        pkpm.skip_field = !pkpmf->evolve;

        app_lw->field_func_ctx = pkpmf->init_ref;
        pkpm.field.init = eval_ic;
        pkpm.field.ctx = &app_lw->field_func_ctx;
      }
    }
  }
  
  app_lw->app = gkyl_pkpm_app_new(&pkpm);
  
  // create Lua userdata ...
  struct pkpm_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct pkpm_app_lw*));
  *l_app_lw = app_lw; // ... point it to the Lua app pointer

  // set metatable
  luaL_getmetatable(L, PKPM_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Apply initial conditions. (time) -> bool
static int
pkpm_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->tstart);
  gkyl_pkpm_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to field. (time) -> bool
static int
pkpm_app_apply_ic_field(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->tstart);
  gkyl_pkpm_app_apply_ic_field(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool
static int
pkpm_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->tstart);
  gkyl_pkpm_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool
static int
pkpm_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_pkpm_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated L2 norm of distribution function (only F_0^2 in 1V). (tm) -> bool
static int
pkpm_app_calc_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_pkpm_app_calc_integrated_L2_f(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool
static int
pkpm_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_pkpm_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool
static int
pkpm_app_write(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_pkpm_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool
static int
pkpm_app_write_field(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_pkpm_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (sidx, time, frame) -> bool
static int
pkpm_app_write_species(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_pkpm_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write diagnostic moments to file (time, frame) -> bool
static int
pkpm_app_write_mom(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_pkpm_app_write_mom(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool
static int
pkpm_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated L2 norm of f to file () -> bool
static int
pkpm_app_write_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_write_integrated_L2_f(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated field energy to file () -> bool
static int
pkpm_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool
static int
pkpm_app_stat_write(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file
static void
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) 
    gkyl_pkpm_app_write(app, tcurr, iot->curr-1);
}

// Run simulation. (num_steps) -> bool. num_steps is optional
static int
pkpm_app_run(lua_State *L)
{
  bool ret_status = true;

  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;
  struct gkyl_pkpm_app *app = app_lw->app;

  double tcurr = app_lw->tstart;
  double tend = app_lw->tend;
  double dt = tend-tcurr;
  long num_steps = luaL_optinteger(L, 2, INT_MAX);

  int nframe = app_lw->nframe;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_pkpm_app_calc_integrated_mom(app, tcurr);
  gkyl_pkpm_app_calc_field_energy(app, tcurr);

  long step = 1;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_pkpm_app_calc_integrated_mom(app, tcurr);
    gkyl_pkpm_app_calc_field_energy(app, tcurr);    
    
    if (!status.success) {
      ret_status = false;
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_pkpm_app_stat_write(app);
  gkyl_pkpm_app_write_integrated_mom(app);
  gkyl_pkpm_app_write_field_energy(app);

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation
static int
pkpm_app_gc(lua_State *L)
{
  struct pkpm_app_lw **l_app_lw = CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor
static struct luaL_Reg pkpm_app_ctor[] = {
  { "new",  pkpm_app_new },
  { 0, 0 }
};

// App methods
static struct luaL_Reg pkpm_app_funcs[] = {
  { "apply_ic", pkpm_app_apply_ic },
  { "apply_ic_field", pkpm_app_apply_ic_field },
  { "apply_ic_species", pkpm_app_apply_ic_species },
  { "calc_integrated_mom", pkpm_app_calc_integrated_mom },
  { "calc_integrated_L2_f", pkpm_app_calc_integrated_L2_f },
  { "calc_field_energy", pkpm_app_calc_field_energy },
  { "write", pkpm_app_write },
  { "write_field", pkpm_app_write_field },
  { "write_species", pkpm_app_write_species },
  { "write_mom", pkpm_app_write_mom },
  { "write_integrated_mom", pkpm_app_write_integrated_mom },
  { "write_integrated_L2_f", pkpm_app_write_integrated_L2_f },
  { "write_field_energy", pkpm_app_write_field_energy },
  { "stat_write", pkpm_app_stat_write },
  { "run", pkpm_app_run },
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App
  do {
    luaL_newmetatable(L, PKPM_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, pkpm_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, pkpm_app_funcs);
    
    luaL_register(L, "G0.PKPM.App", pkpm_app_ctor);
    
  } while (0);

  // Register Species input struct
  do {
    luaL_newmetatable(L, PKPM_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.PKPM.Species", pkpm_species_ctor);
  } while (0);

  // Register Collisions input struct
  do {
    luaL_newmetatable(L, PKPM_COLLISIONS_METATABLE_NM);
    luaL_register(L, "G0.PKPM.Collisions", pkpm_collisions_ctor);
  } while (0);

  // Register Field input struct
  do {
    luaL_newmetatable(L, PKPM_FIELD_METATABLE_NM);
    luaL_register(L, "G0.PKPM.Field", pkpm_field_ctor);
  } while (0);
}

void
gkyl_pkpm_lw_openlibs(lua_State *L)
{
  app_openlibs(L);
}

#endif
