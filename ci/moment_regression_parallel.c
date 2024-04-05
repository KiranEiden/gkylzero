#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

#define RELATIVE_TOLERANCE pow(10.0, -16.0)

int system(const char *command);

void
runTestParallel(const char* test_name, const char* test_name_human, const int test_dimensions, const int test_cuts,
  const int test_output_count, const char test_outputs[][64])
{
  int counter = 0;

  char counter_buffer[128];
  snprintf(counter_buffer, 128, "output_parallel/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
    fclose(counter_ptr);
  }

  counter += 1;

  counter_ptr = fopen(counter_buffer, "w");
  fprintf(counter_ptr, "%d", counter);
  fclose(counter_ptr);

  printf("Running %s...\n", test_name_human);

  char command_buffer1[256];
  snprintf(command_buffer1, 256, "cd ../; rm -rf ./%s-stat.json", test_name);
  system(command_buffer1);
  
  char command_buffer2[256];
  snprintf(command_buffer2, 256, "cd ../; make cuda-build/regression/rt_%s > /dev/null 2>&1", test_name);
  system(command_buffer2);

  char command_buffer3[256];
  if (test_dimensions == 1) {
    snprintf(command_buffer3, 256, "cd ../; mpirun -np %d ./cuda-build/regression/rt_%s -m -M -c %d > ./ci/output_parallel/rt_%s_%d.dat 2>&1",
      test_cuts, test_name, test_cuts, test_name, counter);
  }
  else if (test_dimensions == 2) {
    snprintf(command_buffer3, 256, "cd ../; mpirun -np %d ./cuda-build/regression/rt_%s -m -M -c %d -d %d > ./ci/output_parallel/rt_%s_%d.dat 2>&1",
      test_cuts * test_cuts, test_name, test_cuts, test_cuts, test_name, counter);
  }
  else if (test_dimensions == 3) {
    snprintf(command_buffer3, 256, "cd ../; mpirun -np %d ./cuda-build/regression/rt_%s -m -M -c %d -d %d -e %d > ./ci/output_parallel/rt_%s_%d.dat 2>&1",
      test_cuts * test_cuts * test_cuts, test_name, test_cuts, test_cuts, test_cuts, test_name, counter);
  }
  system(command_buffer3);

  char file_buffer1[128];
  snprintf(file_buffer1, 128, "../%s-stat.json", test_name);
  FILE *file_ptr1 = fopen(file_buffer1, "r");
  if (file_ptr1 == NULL) {
    printf("*** Something catastrophic happened. Test aborting... ***\n");
  }
  else {
    char command_buffer4[256];
    snprintf(command_buffer4, 256, "cd ../; mv ./%s-stat.json ci/output_parallel/%s-stat_%d.json", test_name, test_name, counter);
    system(command_buffer4);
  }

  for (int i = 0; i < test_output_count; i++) {
    char file_buffer2[128];
    snprintf(file_buffer2, 128, "../%s-%s.gkyl", test_name, test_outputs[i]);
    FILE *file_ptr2 = fopen(file_buffer2, "r");
    if (file_ptr2 == NULL) {
      printf("*** Something catastrophic happened. Test aborting... ***\n");
    }
    else {
      char command_buffer5[256];
      snprintf(command_buffer5, 256, "cd ../; mv ./%s-%s.gkyl ci/output/_parallel/%s-%s_%d.gkyl", test_name, test_outputs[i], test_name, test_outputs[i], counter);
      system(command_buffer5);
    }
  }

  printf("Finished %s.\n\n", test_name_human);
}

void
analyzeTestOutputParallel(const char* test_name, const char* test_name_human, const int test_output_count, const char test_outputs[][64])
{
  printf("%s:\n\n", test_name_human);

  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "output_parallel/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
    fclose(counter_ptr);
  }

  int failure[counter + 1];
  int updatecalls[counter + 1];
  int failedsteps[counter + 1];
  double speciesupdate[counter + 1];
  double fieldupdate[counter + 1];
  double sourceupdate[counter + 1];
  double totalupdate[counter + 1];
  int memoryleakcount[counter + 1];
  char *memoryleaks[counter + 1];
  long double averages[counter + 1][test_output_count];

  for (int i = 1; i < counter + 1; i++) {
    char *output;
    long file_size;
    char buffer[128];
    snprintf(buffer, 128, "output_parallel/rt_%s_%d.dat", test_name, i);

    FILE *output_ptr = fopen(buffer, "rb");
    fseek(output_ptr, 0, SEEK_END);
    file_size = ftell(output_ptr);
    rewind(output_ptr);
    output = calloc(file_size, (sizeof(char)));
    fread(output, sizeof(char), file_size, output_ptr);
    fclose(output_ptr);

    updatecalls[i] = 0;
    if (strstr(output, "Number of update calls ") != NULL) {
      char *full_substring = strstr(output, "Number of update calls ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
        substring_index += 1;
      }

      char *end_ptr;
      updatecalls[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    failedsteps[i] = 0;
    if (strstr(output, "Number of failed time-steps ") != NULL) {
      char *full_substring = strstr(output, "Number of failed time-steps ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of failed time-steps ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of failed time-steps ")];
        substring_index += 1;
      }

      char *end_ptr;
      failedsteps[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    speciesupdate[i] = 0.0;
    if (strstr(output, "Species updates took ") != NULL) {
      char *full_substring = strstr(output, "Species updates took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species updates took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species updates took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciesupdate[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    fieldupdate[i] = 0.0;
    if (strstr(output, "Field updates took ") != NULL) {
      char *full_substring = strstr(output, "Field updates took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Field updates took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Field updates took ")];
        substring_index += 1;
      }

      char *end_ptr;
      fieldupdate[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    sourceupdate[i] = 0.0;
    if (strstr(output, "Source updates took ") != NULL) {
      char *full_substring = strstr(output, "Source updates took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Source updates took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Source updates took ")];
        substring_index += 1;
      }

      char *end_ptr;
      sourceupdate[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    totalupdate[i] = 0.0;
    if (strstr(output, "Total updates took ") != NULL) {
      char *full_substring = strstr(output, "Total updates took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
        substring_index += 1;
      }

      char *end_ptr;
      totalupdate[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }
    
    char *temp = output;
    memoryleakcount[i] = 0;
    memoryleaks[i] = (char*)calloc(8192, sizeof(char));
    while (strstr(temp, "0x") != NULL) {
      temp = strstr(temp, "0x");

      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }

      int substring_index = 0;
      int valid_substring = 1;
      while (temp[substring_index] != ' ' && temp[substring_index] != '\n') {
        if (temp[substring_index] != '0' && temp[substring_index] != '1' && temp[substring_index] != '2' && temp[substring_index] != '3' && temp[substring_index] != '4'
          && temp[substring_index] != '5' && temp[substring_index] != '6' && temp[substring_index] != '7' && temp[substring_index] != '8' && temp[substring_index] != '9'
          && temp[substring_index] != 'a' && temp[substring_index] != 'b' && temp[substring_index] != 'c' && temp[substring_index] != 'd' && temp[substring_index] != 'e'
          && temp[substring_index] != 'f' && temp[substring_index] != 'x') {
          valid_substring = 0;
        }

        substring[substring_index] = temp[substring_index];
        substring_index += 1;
      }

      char *temp2 = output;
      int count = 0;
      while (strstr(temp2, substring) != NULL) {
        temp2 = strstr(temp2, substring);

        count += 1;
        temp2 += 1;
      }
      if (count == 1 && valid_substring == 1) {
        memoryleakcount[i] += 1;
        memoryleaks[i] = strcat(memoryleaks[i], substring);
        memoryleaks[i] = strcat(memoryleaks[i], " ");
      }
      
      temp += 1;
    }

    for (int j = 0; j < test_output_count; j++) {
      char *data;
      long data_file_size;
      char data_buffer[256];
      snprintf(data_buffer, 256, "output_parallel/%s-%s_%d.gkyl", test_name, test_outputs[j], i);

      FILE *data_ptr = fopen(data_buffer, "rb");

      if (data_ptr == NULL) {
        failure[i] = 1;
      }
      else {
        fseek(data_ptr, 0, SEEK_END);
        data_file_size = ftell(data_ptr);
        rewind(data_ptr);
        data = calloc(data_file_size, (sizeof(char)));
        fread(data, sizeof(char), data_file_size, data_ptr);
        fclose(data_ptr);

        long long total = 0;
        for (long k = 0; k < data_file_size; k++) {
          total += (long long)abs((int)data[k]);
        }

        averages[i][j] = (long double)total / (long double)data_file_size;
      }
    }
  }

  for (int i = 1; i < counter + 1; i++) {
    printf("Build number: %d\n", i);

    if (failure[i] == 1) {
      printf(ANSI_COLOR_RED "*** Catastrophic test failure ***" ANSI_COLOR_RESET "\n\n");
    }
    else {
      if (i == 1 || failure[i - 1] == 1) {
        printf("Update calls: %d\n", updatecalls[i]);
        printf("Failed time-steps: %d\n", failedsteps[i]);
        printf("Species update time: %f\n", speciesupdate[i]);
        printf("Field update time: %f\n", fieldupdate[i]);
        printf("Source update time: %f\n", sourceupdate[i]);
        printf("Total update time: %f\n", totalupdate[i]);
        if (memoryleakcount[i] != 0) {
          printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", memoryleaks[i]);
        }
        else {
          printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
        }
        printf("Correct: N/A\n\n");
      }
      else {
        if (updatecalls[i] != updatecalls[i - 1]) {
          printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", updatecalls[i]);
        }
        else {
          printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", updatecalls[i]);
        }

        if (failedsteps[i] != failedsteps[i - 1]) {
          printf("Failed time-steps: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", failedsteps[i]);
        }
        else {
          printf("Failed time-steps: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", failedsteps[i]);
        }

        if (speciesupdate[i] > speciesupdate[i - 1]) {
          if (speciesupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Species update time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", speciesupdate[i],
              (((double)speciesupdate[i] / (double)speciesupdate[i - 1]) - 1.0) * 100.0);
          } else {
            printf("Species update time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", speciesupdate[i]);
          }
        }
        else {
          if (speciesupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Species update time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", speciesupdate[i],
              (((double)speciesupdate[i] / (double)speciesupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species update time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", speciesupdate[i]);
          }
        }

        if (fieldupdate[i] > fieldupdate[i - 1]) {
          if (fieldupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Field update time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", fieldupdate[i],
              (((double)fieldupdate[i] / (double)fieldupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Field update time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET, fieldupdate[i]);
          }
        }
        else {
          if (fieldupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Field update time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", fieldupdate[i],
              (((double)fieldupdate[i] / (double)fieldupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Field update time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", fieldupdate[i]);
          }
        }

        if (sourceupdate[i] > sourceupdate[i - 1]) {
          if (sourceupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Source update time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", sourceupdate[i],
              (((double)sourceupdate[i] / (double)sourceupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Source update time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", sourceupdate[i]);
          }
        }
        else {
          if (sourceupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Source update time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", sourceupdate[i],
              (((double)sourceupdate[i] / (double)sourceupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Source update time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", sourceupdate[i]);
          }
        }

        if (totalupdate[i] > totalupdate[i - 1]) {
          if (totalupdate[i  - 1] > pow(10.0, -8.0)) {
            printf("Total update time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", totalupdate[i],
              (((double)totalupdate[i] / (double)totalupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Total update time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", totalupdate[i]);
          }
        }
        else {
          if (totalupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Total update time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", totalupdate[i],
              (((double)totalupdate[i] / (double)totalupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Total update time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", totalupdate[i]);
          }
        }

        if (memoryleakcount[i] != 0) {
          printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", memoryleaks[i]);
        }
        else {
          printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
        }

        int correct = 1;
        for (int j = 0; j < test_output_count; j++) {
          if (fabsl(averages[i][j] - averages[i - 1][j]) > RELATIVE_TOLERANCE) {
            correct = 0;
          }
        }

        if ((updatecalls[i] != updatecalls[i - 1]) || (failedsteps[i] != failedsteps[i - 1]) || (correct != 1)) {
          printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
        }
        else {
          printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
        }
      }
    }
  }
}

void
regenerateTestParallel(const char* test_name, const int test_output_count, const char test_outputs[][64])
{
  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "output_parallel/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
  }
  fclose(counter_ptr);

  for (int i = 1 ; i < counter + 1; i++) {
    char command_buffer[128];
    snprintf(command_buffer, 128, "rm -rf output_parallel/rt_%s_%d.dat", test_name, i);
    system(command_buffer);

    char command_buffer2[128];
    snprintf(command_buffer2, 128, "rm -rf output_parallel/%s-stat_%d.json", test_name, i);
    system(command_buffer2);

    for (int j = 0; j < test_output_count; j++) {
      char command_buffer3[256];
      snprintf(command_buffer3, 256, "rm -rf output_parallel/%s-%s_%d.gkyl", test_name, test_outputs[j], i);
      system(command_buffer3);
    }
  }

  counter_ptr = fopen(counter_buffer, "w");
  fprintf(counter_ptr, "%d", 0);
  fclose(counter_ptr);
}

int
main(int argc, char **argv)
{
  int test_count = 47 - 7;
  char test_names[47 - 7][32] = {
    "10m_burch",
    "10m_burch_grad_closure",
    "10m_gem",
    "10m_gem_grad_closure",
    "10m_lhdi",
    "10m_lhdi_grad_closure",
    "10m_par_firehose",
    "10m_par_firehose_grad_closure",
    "10m_sodshock",
    "10m_sodshock_lax",
    "5m_burch",
    "5m_gem",
    "5m_mom_beach",
    "5m_riem",
    "advect_wv",
    //"advect_wv_mp",
    "burgers_shock",
    //"burgers_shock_mp",
    "coldfluid_beach",
    "coldfluid_clouda",
    "euler_axi_sodshock",
    "euler_axi_vac_riem",
    "euler_bump_in_channel",
    "euler_c2p_sodshock",
    "euler_cart_axi_sodshock",
    "euler_kh_2d",
    "euler_noh_1d",
    "euler_p_perturbation",
    "euler_riem_2d_hll",
    "euler_riem_2d_hllc",
    "euler_riem_2d_lax",
    //"euler_riem_2d_mp",
    "euler_riem_2d_roe",
    "euler_riem_3d",
    "euler_sodshock",
    "euler_sodshock_lax",
    //"euler_sodshock_mp",
    "euler_superwedge",
    "euler_vac",
    "euler_vac_riem_1d",
    //"euler_wave_2d_kep",
    //"euler_wave_2d_mp",
    "euler_wave_2d_wv",
    //"euler_wedge_sodshock",
    "iso_euler_sodshock",
    "iso_euler_sodshock_lax",
    "iso_gem",
  };
  char test_names_human[47 - 7][128] = {
    "Burch et al. Magnetic Reconnection Test (10-moment equations)",
    "Burch et al. Magnetic Reconnection Gradient-Closure Test (10-moment equations)",
    "Geospace Environment Modeling Reconnection Test (10-moment equations)",
    "Geospace Environment Modeling Reconnection Gradient-Closure Test (10-moment equations)",
    "Lower-Hybrid Drift Instability Test (10-moment equations)",
    "Lower-Hybrid Drift Instability Gradient-Closure Test (10-moment equations)",
    "Parallel-Propagating Firehose Instability Test (10-moment equations)",
    "Parallel-Propagating Firehose Instability Gradient-Closure Test (10-moment equations)",
    "Sod-Type Shock Tube Test, with Roe fluxes (10-moment equations)",
    "Sod-Type Shock Tube Test, with Lax fluxes (10-moment equations)",
    "Burch et al. Magnetic Reconnection Test (5-moment equations)",
    "Geospace Environment Modeling Reconnection Test (5-moment equations)",
    "Plasma Wave Beach Test (5-moment equations)",
    "Generalized Brio-Wu Riemann Problem Test (5-moment equations)",
    "Discontinuous Wave Test (linear advection equation)",
    //"Discontinuous Wave Test, with monotonicity-preserving reconstruction (linear advection equation)",
    "Square Wave Test (inviscid Burgers' equation)",
    //"Sinusoidal Wave Test, with monotonicity-preserving reconstruction (inviscid Burgers' equation)",
    "Plasma Wave Beach Test (cold fluid equations)",
    "Dust Cloud Collision Test (cold fluid equations)",
    "2D Sod-Type Shock Tube Test, in axial symmetry/polar coordinates (Euler equations)",
    "2D Sod-Type Vacuum Shock Tube Test, in axial symmetry/polar coordinates (Euler equations)",
    "Supersonic Flow over a Blunt Body Test (Euler equations)",
    "Sod-Type Shock Tube Test, testing computational-to-physical coordinate transformation (Euler equations)",
    "2D Sod-Type Shock Tube Test, in axial symmetry/Cartesian coordinates (Euler equations)",
    "2D Kelvin-Helmholtz Instability Test (Euler equations)",
    "1D Noh Test (Euler equations)",
    "Pressure Perturbation Test (Euler equations)",
    "2D Riemann/quadrant Problem, using Harten-Lax-van Leer solver (Euler equations)",
    "2D Riemann/quadrant Problem, using Harten-Lax-van Leer-with Contact solver (Euler equations)",
    "2D Riemann/quadrant Problem, using Lax-Friedrichs solver (Euler equations)",
    //"2D Riemann/quadrant Problem, with monotonicity-preserving reconstruction (Euler equations)",
    "2D Riemann/quadrant Problem, using Roe solver (Euler equations)",
    "3D Spherical Riemann Problem (Euler equations)",
    "Sod-Type Shock Tube Test, with Roe fluxes (Euler equations)",
    "Sod-Type Shock Tube Test, with Lax fluxes (Euler equations)",
    //"Sod-Type Shock Tube Test, with monotonicity-preserving reconstruction (Euler equations)",
    "2D Superwedge Test (Euler equations)",
    "Sod-Type Vacuum Shock Tube Test (Euler equations)",
    "Double-Rarefaction Riemann Problem Test (Euler equations)",
    //"Smooth Traveling Wave Problem, with kinetic energy-preserving reconstruction (Euler equations)",
    //"Smooth Traveling Wave Problem, with monotonicity-preserving reconstruction (Euler equations)",
    "Smooth Traveling Wave Problem (Euler equations)",
    //"2D Sod-Type Shock Tube Test, with a wedge boundary condition (Euler equations)",
    "Sod-Type Shock Tube Test, with Roe fluxes (isothermal Euler equations)",
    "Sod-Type Shock Tube Test, with Lax fluxes (isothermal Euler equations)",
    "Geospace Environment Modeling Reconnection Test (isothermal Euler equations)",
  };
  int test_dimensions[47 - 7] = { 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, /*1,*/ 1, /*1,*/ 1, 1, 2, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, /*2,*/ 2, 3, 1,
    1, /*1,*/ 2, 1, 1, /*2,*/ /*2,*/ 2, /*2,*/ 1, 1, 2 };
  int test_cuts[47 - 7] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, /*4,*/ 4, /*4,*/ 4, 4, 4, 4, 5, 4, 4, 4, 5, 4, 4, 4, 4, /*4,*/ 4, 1, 4, 4, /*4,*/
    5, 5, 4, /*5,*/ /*5,*/ 5, /*2,*/ 4, 4, 4 };
  int test_output_count[47 - 7] = { 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 4, 4, 3, 4, 1, /*1,*/ 1, /*1,*/ 3, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, /*1,*/ 1, 1,
    1, 1, /*1,*/ 1, 1, 1, /*2,*/ /*1,*/ 1, /*2,*/ 1, 1, 4 };
  char test_outputs[47 - 7][64][64] = {
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "10m_1" },
    { "10m_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "elc_1", "field_1", "ext_em_field_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
    { "q_1" },
    //{ "q_1" },
    { "burgers_1" },
    //{ "burgers_1" },
    { "elc_1", "field_1", "ext_em_field_1" },
    { "cold_1" },
    { "euler_1", "mapc2p" },
    { "euler_1", "mapc2p" },
    { "euler_1", "mapc2p" },
    { "euler_1", "mapc2p" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    //{ "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    //{ "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    //{ "euler_1", "euler-alpha_1" },
    //{ "euler_1" },
    { "euler_1" },
    //{ "euler_1", "mapc2p" },
    { "iso_euler_1" },
    { "iso_euler_1" },
    { "elc_1", "ion_1", "field_1", "ext_em_field_1" },
  };

  system("clear");
  system("mkdir -p output_parallel");

  printf("** Gkeyll Moment App Automated Regression System (Parallel Version) **\n\n");

  if (argc > 1) {
    char *arg_ptr;

    if (strtol(argv[1], &arg_ptr, 10) == 1) {
      for (int i = 0; i < test_count; i++) {
        runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 2) {
      for (int i = 0; i < test_count; i++) {
        analyzeTestOutputParallel(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 3) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          runTestParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_dimensions[strtol(argv[2], &arg_ptr, 10) - 1], test_cuts[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test to run!\n");
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 4) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          analyzeTestOutputParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test results to view!\n");
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 5) {
      for (int i = 0; i < test_count; i++) {
        regenerateTestParallel(test_names[i], test_output_count[i], test_outputs[i]);
        runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 6) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          regenerateTestParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_output_count[strtol(argv[2], &arg_ptr, 10) - 1],
            test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
          runTestParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_dimensions[strtol(argv[2], &arg_ptr, 10) - 1], test_cuts[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test results to (re)generate!\n");
      }
    }
    else {
      printf("Invalid option!\n");
    }
  }
  else {
    while (1) {
      printf("Please select an option to proceed:\n\n");
      printf("1 - Run Full Regression Suite\n");
      printf("2 - View All Regression Results\n");
      printf("3 - Run Specific Regression Test\n");
      printf("4 - View Specific Regression Result\n");
      printf("5 - (Re)generate All Accepted Results\n");
      printf("6 - (Re)generate Specific Accepted Result\n");
      printf("7 - Exit\n");

      int option;
      scanf("%d", &option);
      printf("\n");

      if (option == 1) {
        for (int i = 0; i < test_count; i++) {
          runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 2) {
        for (int i = 0; i < test_count; i++) {
          analyzeTestOutputParallel(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 3) {
        printf("Please select the test you wish to run:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          runTestParallel(test_names[option2 - 1], test_names_human[option2 - 1], test_dimensions[option2 - 1], test_cuts[option2 - 1],
            test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 4) {
        printf("Please select the test whose results you wish to view:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          analyzeTestOutputParallel(test_names[option2 - 1], test_names_human[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 5) {
        for (int i = 0; i < test_count; i++) {
          regenerateTestParallel(test_names[i], test_output_count[i], test_outputs[i]);
          runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 6) {
        printf("Please select the test whose accepted result you wish to (re)generate:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          regenerateTestParallel(test_names[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
          runTestParallel(test_names[option2 - 1], test_names_human[option2 - 1], test_dimensions[option2 - 1], test_cuts[option2 - 1],
            test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 7) {
        break;
      }
      else {
        printf("Invalid selection!\n\n");
      }
    }
  }

  return 0;
}