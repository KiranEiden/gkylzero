// Test the FEM Poisson solver with a nonzero kSq,
// essentially solving the Helmholtz equation.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson.h>

void evalFunc2x_dirichletx_periodicy_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double m = 5.;
  double kSq = 0.3;
  // RHSsource:
  //   ( 4*gxx*M_PI - 90*gxx*x - kSq*(4*M_PI - 5*x)*(M_PI - 3*x)*(M_PI + x)
  //    +10*gxy*m*(3*pow(M_PI,2) + 2*M_PI*x - 9*pow(x,2))*cos(m*y)
  //    -5*(-2*gxx*(M_PI - 9*x) + (kSq - gyy*pow(m,2))*(M_PI - 3*x)*(M_PI - x)*(M_PI + x))*sin(m*y))/15.;
  fout[0] = ((x-4.*M_PI/5.) + (x-M_PI)*sin(m*y))*(x+M_PI)*(x-M_PI/3.);
}
void evalFunc2x_dirichletx_periodicy(double t, const double *xn, double* restrict fout, void *ctx)
{
  // These values have to match those in the test below.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double x = xn[0], y = xn[1];
  double m = 5.;
  double kSq = 0.3;
  // Expected solution: ((x-4.*M_PI/5.) + (x-M_PI)*sin(m*y))*(x+M_PI)*(x-M_PI/3.);
  fout[0] = 
    ( 4*gxx*M_PI - 90*gxx*x - kSq*(4*M_PI - 5*x)*(M_PI - 3*x)*(M_PI + x)
     +10*gxy*m*(3*pow(M_PI,2) + 2*M_PI*x - 9*pow(x,2))*cos(m*y)
     -5*(-2*gxx*(M_PI - 9*x) + (kSq - gyy*pow(m,2))*(M_PI - 3*x)*(M_PI - x)*(M_PI + x))*sin(m*y) )/15.;
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

// Apply periodic BCs in one direction.
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, sgr.lower_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(buff->data, fld, sgr.upper_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.lower_ghost[dir]);
}

void
test_2x(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  // Determinant of g tensor has to be >0. Diagonal entries have to be  >0.
  double gxx = 3.0;
  double gxy = 2.44948;
  double gyy = 2.0;

  double kSq = 0.3;

  double lower[] = {-M_PI,-M_PI}, upper[] = {M_PI,M_PI};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // Local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // Skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob, *projob_sol;
  if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
             (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_periodicy, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc2x_dirichletx_periodicy_sol, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = dim+ceil((pow(3.,dim-1)-dim)/2);
  struct gkyl_array *eps = use_gpu? mkarr_cu(epsnum*basis.num_basis, localRange_ext.volume)
                                  : mkarr(epsnum*basis.num_basis, localRange_ext.volume);
  struct gkyl_array *kSqFld;  // kSq field multiplying phi in Helmholtz equation.
  // Device copies:
  struct gkyl_array *rho_cu, *phi_cu;
  if (use_gpu) {
    rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    kSqFld = mkarr_cu(basis.num_basis, localRange_ext.volume);
  } else {
    kSqFld = mkarr(basis.num_basis, localRange_ext.volume);
  }

  // Project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  if (use_gpu) gkyl_array_copy(rho_cu, rho);
//  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_poisson_vareps_2x_rho_1.gkyl");

  // Project the expected solution.
  if (projob_sol != NULL) {
    struct gkyl_array *phisol = mkarr(basis.num_basis, localRange_ext.volume);
    gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol);
//    gkyl_grid_sub_array_write(&grid, &localRange, phisol, "ctest_fem_helmholtz_2x_phisol.gkyl");
    gkyl_array_release(phisol);
    gkyl_proj_on_basis_release(projob_sol);
  }

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, gxx*dg0norm, 0*basis.num_basis);
  gkyl_array_shiftc(eps, gxy*dg0norm, 1*basis.num_basis);
  gkyl_array_shiftc(eps, gyy*dg0norm, 2*basis.num_basis);
//  gkyl_grid_sub_array_write(&grid, &localRange, eps, "ctest_fem_poisson_vareps_2x_eps_1.gkyl");

  // Project the kSq onto the basis.
  gkyl_array_shiftc(kSqFld, kSq*dg0norm, 0);

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&grid, basis, &bcs, 0, eps, kSqFld, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_poisson_set_rhs(poisson, rho_cu);
  else
    gkyl_fem_poisson_set_rhs(poisson, rho);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_poisson_solve(poisson, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_poisson_solve(poisson, phi);
  }
//  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_helmholtz_2x_phi_8x8_p2.gkyl");

  if (poly_order == 1) {
    if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8:
      const double sol[256] = {
        2.2723820361020547e+01, 1.3119603802451259e+01, -4.0251978322074500e+00, -2.3239490519656978e+00, 3.9871828304012533e+00, 2.3020010804407463e+00,
        -6.7924048898040432e+00, -3.9215967915735011e+00, 1.5832165717984553e+01, 9.1407051391335976e+00, 1.3631108948496774e+01, 7.8699244207680055e+00,
        1.7817467802904730e+01, 1.0286919832284779e+01, -1.2484894255345385e+01, -7.2081570591280668e+00, 3.1648437814200294e+00, 1.8272234091458350e+00,
        4.0251978322068140e+00, 2.3239490519662831e+00, 2.1901481312039152e+01, 1.2644826131156655e+01, 6.7924048898045832e+00, 3.9215967915730929e+00,
        1.0056498424454976e+01, 5.8061220724639435e+00, -1.3631108948497818e+01, -7.8699244207675161e+00, 8.0711963395349660e+00, 4.6599073793124566e+00,
        1.2484894255346525e+01, 7.2081570591273998e+00, 5.2068465407884489e+01, 3.8225349146296583e+00, 1.5595761024353239e+00, 5.5483197864949103e+00,
        1.4583976837804151e+01, 3.8160607922803313e+00, -2.3201255668806660e+01, -5.5520576227760605e+00, 2.8528195952778681e+01, -1.8106486654355789e+00,
        3.1251954328476330e+01, 2.3034754027119662e+00, 4.6292580733758555e+01, 6.1531942835486744e+00, -2.0995681993189507e+01, 2.2944514676679968e+00,
        7.2257277343075419e+00, 5.1732903420160214e-01, -1.5595761024350800e+00, -5.5483197864949876e+00, 4.4710216304388261e+01, 5.2380315655093990e-01,
        2.3201255668806635e+01, 5.5520576227761431e+00, 3.0765997189413419e+01, 6.1505126142670354e+00, -3.1251954328476486e+01, -2.3034754027119408e+00,
        1.3001612408433161e+01, -1.8133303347172267e+00, 2.0995681993189443e+01, -2.2944514676680265e+00, 5.2440522267967133e+01, -3.6077277863070969e+00,
        1.5873676018183975e+01, 2.7159296530697716e+00, 2.3672467455535674e+01, 1.4311817123943515e+00, -3.2482920208203183e+01, 1.9328610242767283e-01,
        1.9482965775977302e+01, -3.4116174120228222e+00, 3.0064110285739766e+01, -2.9892774805412126e+00, 5.4175870683224048e+01, -1.6017747091249026e+00,
        -1.0034152298570771e+01, 4.0341906522502073e+00, 9.3021957248332363e+00, 6.8152031909205002e-01, -1.5873676018183936e+01, -2.7159296530698160e+00,
        3.8070250537264926e+01, -4.3573891796094975e+00, 3.2482920208203282e+01, -1.9328610242768393e-01, 4.2259752216823379e+01, 4.8540994480772282e-01,
        -3.0064110285739815e+01, 2.9892774805412503e+00, 7.5668473095764091e+00, -1.3244327580901019e+00, 1.0034152298570683e+01, -4.0341906522501887e+00,
        3.5044168186816982e+01, -6.4360619253629894e+00, 2.0348561428026951e+01, -1.3235335710419840e-01, 2.3784564026096898e+01, -1.3664627272055949e+00,
        -2.6849296921187321e+01, 3.0592878188438819e+00, 7.8025527662067518e+00, -3.3320721834142146e+00, 1.7622078418098109e+01, -4.1941329675076027e+00,
        4.1664134004662145e+01, -5.6218798302725270e+00, 1.9279146251104966e+00, 2.8721119062014928e+00, 9.7586378337295834e+00, -4.1799334465126453e-01,
        -2.0348561428026986e+01, 1.3235335710419946e-01, 2.1018241994449660e+01, -5.4875925428086951e+00, 2.6849296921187356e+01, -3.0592878188439032e+00,
        3.7000253254339917e+01, -3.5219830866001067e+00, -1.7622078418098084e+01, 4.1941329675076053e+00, 3.1386720158845236e+00, -1.2321754397417586e+00,
        -1.9279146251105217e+00, -2.8721119062014746e+00, 1.2429370149904974e+01, -6.6205978089171778e+00, 1.6959032197075697e+01, -1.8245922568116115e+00,
        1.6916080202993015e+01, -2.5990582569880232e+00, -1.4368628940149517e+01, 4.1464291996746363e+00, -2.1492113622285096e+00, -2.4135815150495454e+00,
        3.3612777227901676e+00, -4.0393441527880407e+00, 2.0326472486816517e+01, -6.6974247909327680e+00, 9.6150643976769050e+00, 1.5660660842906933e+00,
        7.6063471291079514e+00, -8.2463227303636466e-01, -1.6959032197075729e+01, 1.8245922568116157e+00, 3.1196370760198544e+00, -4.8461718249655137e+00,
        1.4368628940149517e+01, -4.1464291996746381e+00, 2.2184928641241402e+01, -5.0316485669040114e+00, -3.3612777227901547e+00, 4.0393441527880318e+00,
        -2.9075520780356889e-01, -7.4780529102079274e-01, -9.6150643976768855e+00, -1.5660660842906868e+00, -6.2750640026612921e+00, -4.1784122841065878e+00,
        8.2564814069615000e+00, -3.1998277844972116e+00, 7.7716710026903790e+00, -2.6804688563867942e+00, -1.4659517038637268e-01, 4.0646658257220656e+00,
        -6.4238012029651843e+00, -5.4354080175013225e-02, -8.0491645288227200e+00, -2.5484777527533398e+00, -3.9510689078868699e-01, -5.2661846407617965e+00,
        1.1529832812819787e+01, -4.6057402437215789e-01, 5.2745040551848454e+00, -5.2165795340110488e-01, -8.2564814069615178e+00, 3.1998277844972169e+00,
        -8.7722309501668576e+00, -2.0196013811208897e+00, 1.4659517038636960e-01, -4.0646658257220656e+00, 5.4232412554886942e+00, -4.6457161573326706e+00,
        8.0491645288227147e+00, 2.5484777527533398e+00, -6.0545305668776983e-01, 5.6611440325410312e-01, -1.1529832812819762e+01, 4.6057402437215300e-01,
        -1.0802498418910655e+01, 1.5644968051466410e+00, -2.4972452249550914e+00, -3.0088391812314819e+00, 6.1266477647347584e-01, -1.4527853154497430e+00,
        9.0877927686536299e+00, 1.2668105366841138e+00, -1.5819159823365418e+00, 2.8498178156902303e+00, -1.0354834560510966e+01, 1.2172985392956770e+00,
        -9.8934733048884986e+00, -2.1769976528043039e-01, 5.5561547029505878e+00, -2.9883306404129977e+00, 4.0554245439156702e+00, -1.8217793059371093e-01,
        2.4972452249550905e+00, 3.0088391812314845e+00, -7.3597386514684677e+00, 2.8351041900026783e+00, -9.0877927686536317e+00, -1.2668105366841140e+00,
        -5.1651578926584536e+00, -1.4674989411372898e+00, 1.0354834560510968e+01, -1.2172985392956739e+00, 3.1463994298935098e+00, 1.6000186398333658e+00,
        -5.5561547029505860e+00, 2.9883306404129919e+00, 4.2219852160212206e+00, 7.1098928659164056e+00, -3.8543537793259874e+00, 2.2253121920459011e+00,
        7.3165238468922222e+00, 5.3232601543653217e+00, 5.6409864908770562e+00, -3.2568250690029115e+00, 9.9453970814568979e+00, 3.8054794847240023e+00,
        -4.1232058212357376e+00, 2.3805339908146852e+00, 3.1330702685049392e+00, 7.7385782039895217e+00, 1.9010710197021807e-01, -1.0975838649736505e-01,
        1.0138282004134755e+01, 3.6941173226933226e+00, 3.8543537793259817e+00, -2.2253121920459074e+00, 7.0437433732637391e+00, 5.4807500342443989e+00,
        -5.6409864908770571e+00, 3.2568250690029119e+00, 4.4148701386990661e+00, 6.9985307038857156e+00, 4.1232058212357412e+00, -2.3805339908146874e+00,
        1.1227196951651038e+01, 3.0654319846202047e+00, -1.9010710197021397e-01, 1.0975838649737209e-01,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
      TEST_MSG("This BC combination is not available");
    }
  } else if (poly_order == 2) {
    if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1]==GKYL_POISSON_PERIODIC && bcs.up_type[1]==GKYL_POISSON_PERIODIC)) {
      // Solution with N=8x8:
      const double sol[512] = {
        7.6005338241730174e+00, 2.5360186868119099e+00, 5.2094857775142778e+00, 2.5292000192143851e+00, -1.4346704312888237e+00, 6.1418092466141374e+00,
        -3.7064295453898721e-01, 3.5459752218440470e+00, 2.0875531459905702e+01, 1.1363293721363824e+01, -1.0125609749462622e+01, -4.6960951122899086e+00,
        -5.3385201397583382e-01, -3.6029634305647162e+00, 8.9073071026561490e-01, -2.0801719065169442e+00, 1.3400405151192054e+01, 7.0513069680342060e+00,
        9.1102888574731029e+00, 4.1120813787797399e+00, -5.3092447758244410e-01, -1.0464494983749311e+00, -8.8904049634082838e-01, -6.0416789958047501e-01,
        1.0696832521693558e+01, 4.3221020808140258e+00, -2.7582843099125176e+00, -1.1192661431623625e+00, -1.4358830565672702e+00, 5.0828665035055511e+00,
        3.6656241715847870e-01, 2.9345943440537043e+00, 2.1995387909904061e+01, 1.2493767400145195e+01, -5.2094857775141001e+00, -2.5292000192145050e+00,
        -1.5900589717447083e-01, -6.1418092466140699e+00, 3.7064295453903850e-01, -3.5459752218440785e+00, 8.7203902741715691e+00, 3.6664923655932471e+00,
        1.0125609749462779e+01, 4.6960951122898367e+00, -1.0598243144874206e+00, 3.6029634305649574e+00, -8.9073071026564210e-01, 2.0801719065167799e+00,
        1.6195516582884629e+01, 7.9784791189230644e+00, -9.1102888574733676e+00, -4.1120813787798260e+00, -1.0627518508807672e+00, 1.0464494983753512e+00,
        8.8904049634088089e-01, 6.0416789958009898e-01, 1.8899089212383547e+01, 1.0707684006142872e+01, 2.7582843099124466e+00, 1.1192661431626412e+00,
        -1.5779327189598261e-01, -5.0828665035055813e+00, -3.6656241715855459e-01, -2.9345943440537652e+00, 1.0198582668159245e+01, 4.3022079578716810e-01,
        8.1683343508162949e+00, -5.3199983739823375e-01, -2.9891397269858605e-01, 1.3834673774896498e+01, -1.4685565490202868e-01, 8.9550218439909635e-01,
        4.3559800646084057e+01, 1.3207369343973014e+00, -1.5656805206259619e+01, 1.2891233386180193e-01, -8.5355769950612637e-01, -6.8854610496074171e+00,
        -1.7336805007464207e-01, 1.8502102254833480e-01, 3.4409804580016520e+01, 5.1509121431144482e+00, 1.3973731915309614e+01, 3.4969026649379437e-01,
        -4.7481525459127716e-01, -4.0971613753500629e+00, 3.9203510259976782e-01, -1.1571614238111050e+00, 1.3988635134384170e+01, -1.1562897219286308e+00,
        -4.1050359853373699e+00, -6.2344905136715567e-01, -4.5579423002864322e-01, 1.2679722233858552e+01, -3.8105330894826900e-01, 1.4514523568603499e+00,
        5.2018525989982813e+01, 3.9332654875229034e+00, -8.1683343508162256e+00, 5.3199983739822099e-01, -8.6143646585017919e-01, -1.3834673774896491e+01,
        1.4685565490202382e-01, -8.9550218439910123e-01, 1.8657308012058241e+01, 3.0427493489127655e+00, 1.5656805206259653e+01, -1.2891233386176920e-01,
        -3.0679273904265281e-01, 6.8854610496073958e+00, 1.7336805007463907e-01, -1.8502102254832178e-01, 2.7807304078125789e+01, -7.8742585980429280e-01,
        -1.3973731915309783e+01, -3.4969026649371976e-01, -6.8553518395750324e-01, 4.0971613753499270e+00, -3.9203510259976676e-01, 1.1571614238111605e+00,
        4.8228473523757707e+01, 5.5197760052387306e+00, 4.1050359853374347e+00, 6.2344905136706252e-01, -7.0455620852012613e-01, -1.2679722233858572e+01,
        3.8105330894827566e-01, -1.4514523568602842e+00, 1.2586256325308968e+01, 8.4107295768330115e-01, 5.1957341434021052e+00, -5.8594003257720051e-01,
        -3.8197418672731953e-01, 1.3691090896679075e+01, 3.1657907214436448e-01, -9.7839979778894859e-01, 3.8023826745903961e+01, -3.6600968274003298e+00,
        -1.4859549348382227e+01, 2.8969588802972629e-01, -1.8985528828992038e-01, -5.1704955018413994e+00, -2.0565917441511680e-01, 8.0511479810530950e-01,
        4.2387458583333789e+01, -6.9453852908623692e-01, 1.5818842075832473e+01, 1.7624817876190521e-01, -5.9063842674754496e-01, -6.3789060337859462e+00,
        -2.5733078460055018e-02, -1.6020446895884380e-01, 1.0778780837042301e+01, -3.8730150948623548e-01, -7.5116716562984580e+00, -5.3894845277832537e-01,
        -2.1596437520770667e-01, 1.4191630927925010e+01, 2.4205124297494790e-01, -5.7855146535094415e-01, 5.1116569757097245e+01, -3.7873585677816437e+00,
        -5.1957341434020705e+00, 5.8594003257718807e-01, -3.4505036190696492e-01, -1.3691090896679071e+01, -3.1657907214437242e-01, 9.7839979778894970e-01,
        2.5678999336502432e+01, 7.1381121730195662e-01, 1.4859549348382302e+01, -2.8969588802972868e-01, -5.3716926034437829e-01, 5.1704955018414065e+00,
        2.0565917441511788e-01, -8.0511479810530684e-01, 2.1315367499072842e+01, -2.2517470810120970e+00, -1.5818842075832455e+01, -1.7624817876187754e-01,
        -1.3638612188675606e-01, 6.3789060337859178e+00, 2.5733078460052457e-02, 1.6020446895884879e-01, 5.2924045245364006e+01, -2.5589841006119829e+00,
        7.5116716562983319e+00, 5.3894845277831338e-01, -5.1106017342658194e-01, -1.4191630927924978e+01, -2.4205124297493713e-01, 5.7855146535090951e-01,
        1.1545061003011186e+01, -1.3738720815870096e+00, 4.3612849575033330e+00, -3.4498508344333356e-01, -3.2904190977538550e-01, 9.4721679538244210e+00,
        -3.1335401350579467e-02, -1.4573964989584725e+00, 2.4423142752336549e+01, -3.9909454407253895e+00, -1.2826790767726827e+01, 1.4098501079087633e+00,
        -3.3916355383896001e-02, -2.8444378924773055e+00, 2.0172621962324291e-01, 5.3783518881162773e-01, 3.1896767280843616e+01, -4.7603125051136201e+00,
        1.3778536507938240e+01, -1.6488440600743914e+00, -1.2436809238044537e-01, -5.4495253089551019e+00, -2.5394855432687014e-01, 6.9678268061958215e-01,
        8.4493843632199024e+00, -1.0551898090740577e+00, -6.6590024314522509e+00, 9.2196752408683447e-01, -2.9157557357122887e-01, 1.0551230492896986e+01,
        1.5741127005085176e-01, -1.5232347057705165e+00, 3.4135366758970932e+01, -5.5256575116810733e+00, -4.3612849575033952e+00, 3.4498508344331402e-01,
        3.5343251055574419e-02, -9.4721679538244246e+00, 3.1335401350589494e-02, 1.4573964989584691e+00, 2.1257285009645475e+01, -2.9085841525427933e+00,
        1.2826790767726866e+01, -1.4098501079087915e+00, -2.5978230333590824e-01, 2.8444378924773228e+00, -2.0172621962325094e-01, -5.3783518881162362e-01,
        1.3783660481138687e+01, -2.1392170881545916e+00, -1.3778536507938163e+01, 1.6488440600744003e+00, -1.6933056633937377e-01, 5.4495253089550832e+00,
        2.5394855432686930e-01, -6.9678268061957960e-01, 3.7231043398762452e+01, -5.8443397841940703e+00, 6.6590024314521994e+00, -9.2196752408679683e-01,
        -2.1230851485945852e-03, -1.0551230492897009e+01, -1.5741127005085417e-01, 1.5232347057705182e+00, 5.2910994659911337e+00, -1.6037930706593944e+00,
        1.9472396632635633e+00, -1.0281253693285326e+00, 1.6132520626729868e-01, 4.0722279539752053e+00, -1.5348308831611070e-02, -1.6602603135622993e+00,
        9.0527570817568144e+00, -4.8428263726504470e+00, -5.8202254246033549e+00, 2.2725363186585170e+00, -2.6811389031212618e-03, -6.1003834416465774e-01,
        -7.9340798701337728e-02, 7.5219599188385800e-01, 1.5668276072738051e+01, -4.3443080031298349e+00, 6.2838020682795177e+00, -2.1857263135037730e+00,
        8.0825089942821546e-02, -3.2095034540899299e+00, 1.2755314240455373e-01, 5.9649454027746884e-01, 2.5508617777908702e+00, -1.8102861404069106e+00,
        -3.0664126836254448e+00, 8.1854747753430213e-01, 1.2673579373666191e-01, 5.1489616574218964e+00, -1.0104658521049491e-01, -1.5957666606257512e+00,
        1.4486167985990903e+01, -5.8924525955400213e+00, -1.9472396632636046e+00, 1.0281253693285521e+00, -2.1697975072643833e-02, -4.0722279539752160e+00,
        1.5348308831611505e-02, 1.6602603135622993e+00, 1.0724510370225014e+01, -2.6534192935489336e+00, 5.8202254246033025e+00, -2.2725363186585255e+00,
        1.4230837009778921e-01, 6.1003834416466907e-01, 7.9340798701344584e-02, -7.5219599188386543e-01, 4.1089913792438120e+00, -3.1519376630696252e+00,
        -6.2838020682794600e+00, 2.1857263135037539e+00, 5.8802141251856774e-02, 3.2095034540899290e+00, -1.2755314240455398e-01, -5.9649454027746218e-01,
        1.7226405674191181e+01, -5.6859595257925646e+00, 3.0664126836254813e+00, -8.1854747753429369e-01, 1.2891437458003635e-02, -5.1489616574219088e+00,
        1.0104658521048825e-01, 1.5957666606257566e+00, 2.1875877255895033e+00, -6.9629386298533269e-02, -1.1122354167860504e+00, -4.2477508077066245e-01,
        2.5303038897513441e-01, -1.3331215683293596e+00, 2.2747869267361009e-01, -1.4605196882042557e+00, -5.3879509883823395e+00, -2.9417048701490547e+00,
        3.9588849814982607e-01, 1.2536655577031333e+00, 4.2552786643030904e-01, 1.4998824005980973e+00, -1.2788722763033364e-01, 4.6596731807369385e-01,
        1.3704046061212702e+00, -3.8557699977809943e+00, 5.5236453351529080e-01, -1.3481757536131740e+00, 1.2327448949887937e-01, -7.8803246456118026e-01,
        -4.6618840904500726e-02, 8.0154238736182382e-01, -6.1181482099308648e-01, 3.0898878645883948e-01, -1.1770499128207832e+00, 6.5294287751918489e-01,
        3.7822783697318579e-01, -3.8543620162538889e-01, 1.9381622469958776e-01, -1.5995194330976867e+00, -4.5668888180714626e+00, -4.6668044425936221e+00,
        1.1122354167860498e+00, 4.2477508077066312e-01, 3.1992273213401451e-01, 1.3331215683293527e+00, -2.2747869267361268e-01, 1.4605196882042586e+00,
        3.0086498959003167e+00, -1.7947289587430753e+00, -3.9588849814985644e-01, -1.2536655577031235e+00, 1.4742525467883247e-01, -1.4998824005980997e+00,
        1.2788722763033192e-01, -4.6596731807369429e-01, -3.7497056986033579e+00, -8.8066383111112301e-01, -5.5236453351528647e-01, 1.3481757536131691e+00,
        4.4967863161026617e-01, 7.8803246456118659e-01, 4.6618840904504730e-02, -8.0154238736182637e-01, -1.7674862714889301e+00, -5.0454226153509847e+00,
        1.1770499128208101e+00, -6.5294287751919011e-01, 1.9472528413596699e-01, 3.8543620162538739e-01, -1.9381622469958751e-01, 1.5995194330976885e+00,
        2.1366515781456719e+00, 1.9064759605045620e-01, -8.2574641291077466e-01, 1.0988937788293476e-01, 3.6955002880581328e-01, -4.1196685239664932e+00,
        -1.4455250146925050e-01, -1.4829394674237523e-01, -7.0526350785989278e+00, 2.1080735609429442e+00, 3.0422153923089352e+00, 7.4123101757044585e-01,
        5.2426918830117542e-01, 1.7067748907949338e+00, 2.3387964985269036e-01, -3.4651788316523829e-01, -7.8031021871712305e+00, -8.1626706528285342e-01,
        -3.4765958545526261e+00, -1.1581483357827167e+00, 6.0684709824848970e-01, 1.7059243254864720e+00, -1.8620327131549172e-01, 6.3834423671948293e-01,
        2.4475052326316993e+00, 1.4019491444316479e+00, 1.8744336160896424e+00, 8.9663806613325436e-01, 3.3534513855320786e-01, -4.1193162082800505e+00,
        2.9451541799910250e-02, -5.5623719386614334e-01, -1.1298575694442931e+01, 1.1892583226032802e+00, 8.2574641291077422e-01, -1.0988937788292615e-01,
        6.3672898221781615e-01, 4.1196685239664861e+00, 1.4455250146925511e-01, 1.4829394674237223e-01, -2.1092890376983897e+00, -7.2816764228921660e-01,
        -3.0422153923089574e+00, -7.4123101757045151e-01, 4.8200982272245790e-01, -1.7067748907949318e+00, -2.3387964985269274e-01, 3.4651788316524151e-01,
        -1.3588219291261323e+00, 2.1961729839365658e+00, 3.4765958545526328e+00, 1.1581483357827163e+00, 3.9943191277513734e-01, -1.7059243254864613e+00,
        1.8620327131549039e-01, -6.3834423671947771e-01, -1.1609429348929011e+01, -2.2043225777933582e-02, -1.8744336160896267e+00, -8.9663806613325703e-01,
        6.7093387247041492e-01, 4.1193162082800558e+00, -2.9451541799911076e-02, 5.5623719386614556e-01, 7.4732049141385701e+00, 3.8230630856625507e+00,
        -7.0106178850855883e-01, 2.7673600642057300e-01, 1.0919777491929925e+00, -2.1882605870896055e+00, 9.9165574273978152e-02, 1.2633928390132259e+00,
        6.0202620362782389e+00, 5.4171978776716267e+00, 3.4015367807192884e+00, -1.3997961114424651e+00, 5.0694205744484344e-01, 5.5329415571072216e-01,
        -4.3693608839031134e-01, -3.1944452974067922e-01, 2.8879154137810086e+00, 7.0427055726954677e+00, -4.1094376596956188e+00, 1.7028746389384624e+00,
        6.4865720645202729e-01, 1.4057844881016739e+00, 5.1875536781785714e-01, -8.1163005262812660e-01, 8.7706653672304427e+00, 3.1497557526417608e+00,
        2.4100856913490154e+00, -1.0084322979655316e+00, 1.0332774124805153e+00, -2.5413736445578960e+00, -2.9669478833153989e-01, 1.4672627577968969e+00,
        3.5835472215099453e+00, 7.0297104907755887e+00, 7.0106178850862011e-01, -2.7673600642055746e-01, 3.4762715174511066e-01, 2.1882605870896032e+00,
        -9.9165574273982746e-02, -1.2633928390132199e+00, 5.0364900993703268e+00, 5.4355756987665860e+00, -3.4015367807193400e+00, 1.3997961114424662e+00,
        9.3266284349326589e-01, -5.5329415571073826e-01, 4.3693608839031900e-01, 3.1944452974066434e-01, 8.1688367218673896e+00, 3.8100680037427059e+00,
        4.1094376596955895e+00, -1.7028746389384872e+00, 7.9094769448608682e-01, -1.4057844881016783e+00, -5.1875536781786191e-01, 8.1163005262811183e-01,
        2.2860867684179671e+00, 7.7030178237963876e+00, -2.4100856913489959e+00, 1.0084322979655396e+00, 4.0632748845759326e-01, 2.5413736445578783e+00,
        2.9669478833154184e-01, -1.4672627577969128e+00,
      };
      for (int j=0; j<cells[0]; j++) {
        for (int k=0; k<cells[1]; k++) {
          int idx0[] = {j+1,k+1};
          long linidx = gkyl_range_idx(&localRange, idx0);
          const double *phi_p = gkyl_array_cfetch(phi, linidx);
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[(j*cells[1]+k)*basis.num_basis+m], phi_p[m], 1e-10) );
            TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[(j*cells[1]+k)*basis.num_basis+m], j, k);
            TEST_MSG("Produced: %.13e", phi_p[m]);
          }
        }
      }
    } else {
      TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
      TEST_MSG("This BC combination is not available");
    }
  } else {
    TEST_CHECK( gkyl_compare(1., 2., 1e-10) );
    TEST_MSG("This poly_order is not available");
  }

  gkyl_fem_poisson_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(phi);
  gkyl_array_release(kSqFld);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }
}

void test_2x_p1_dirichletx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(1, &cells[0], bc_tv, false);
}

void test_2x_p2_dirichletx_periodicy() {
  int cells[] = {8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[1] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 4.*pow(M_PI,3)/15.;
  test_2x(2, &cells[0], bc_tv, false);
}

TEST_LIST = {
  // 2x tests
  { "test_2x_p1_dirichletx_periodicy", test_2x_p1_dirichletx_periodicy },
  { "test_2x_p2_dirichletx_periodicy", test_2x_p2_dirichletx_periodicy },
#ifdef GKYL_HAVE_CUDA
#endif
  { NULL, NULL },
};
