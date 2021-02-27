// README
// 
// Gray-Scott Reaction-Diffusion
// An implementation in Processing by Kaylah Facey.
//
// Commands:
// i - Initialize the system with a fixed rectangular region that has specific u and v concentrations.
// space bar - Start or stop the simulation (toggle between these).
// u - At each timestep, draw values for u at each cell (default).
// v - At each timestep, draw values for v at each cell (rather than u).
// d - Toggle between performing diffusion alone or reaction-diffusion (reaction-diffusion is default).
// p - Toggle between constant f, k for each cell (default) and spatially-varying f, k.
// 1 - Set parameters for spots (k = 0.0625, f = 0.035)
// 2 - Set parameters for stripes (k = 0.06, f = 0.035)
// 3 - Set parameters for spiral waves (k = 0.0475, f = 0.0118)
// 4 - Set parameters for a custom pattern (k = , f = )
// (mouse click) - Print values for u and v at cell.  If in spatially-varying parameter mode, also print values for k and f at the cell.
//
// Custom Commands:
// t - Set boundary condition to toroidal.
// z - Set boundary condition to zero-derivative.
// f - Set boundary condition to fixed.


// Simulation Parameters.
float ru = 0.082;
float rv = 0.041;
float f = 0;
float k = 0;

float[] f_case = new float[]{.035, .035, .0118, 0};
float[] k_case = new float[]{.0625, .06, .0475, 0};

// Colours.
color red = color(255, 0, 0);

// Drawing
int cells_per_side = 100;                        // cells per side of grid.
int cell_size = 8;                               // size of a cell.
boolean running = false;                         // Whether the simulation is running (true) or stopped (false).
boolean draw_u = true;                           // Whether to draw values for u at each cell (true) 
                                                 // or to draw values for v at each cell instead (false).
boolean diffusion_only = false;                  // Whether to perform reaction-diffusion (false) or diffusion only (true).
boolean constant_fk = true;                      // Toggle between constant f, k for each cell (true) and spatially-varying f, k (false).

// Boundary condition values.
int toroidal = 0;
int zero_derivative = 1;
int fixed = 2;
int boundary_condition = toroidal;               

// Cell state.
int curr_case = 0;                                                       // The current case - 1.
float[][] u_curr = new float[cells_per_side][cells_per_side];            // Current concentrations u
float[][] u_next = new float[cells_per_side][cells_per_side];            // Future concentrations u
float[][] v_curr = new float[cells_per_side][cells_per_side];            // Current concentrations v
float[][] v_next = new float[cells_per_side][cells_per_side];            // Future concentrations v
float[][] k_curr = new float[cells_per_side][cells_per_side];            // k values
float[][] f_curr = new float[cells_per_side][cells_per_side];            // f values

// Fixed region of varying u, v
int fixed_x = cells_per_side/2 - 5;
int fixed_y = cells_per_side/2 - 5;

// Track min and max values of u and v.
float min_u = 1;
float max_u = 0;
float min_v = 1;
float max_v = 0;

void setup() {
  size(800, 800);  // cell_size * cells_per_side
  set_case(0);
}

// Update parameters for a new case.
void set_case(int new_case) {
  println("Set case to " + (new_case + 1) + ".");
  curr_case = new_case;
  f = f_case[curr_case];
  k = k_case[curr_case];
  initialize();
}

// Fill grid based on current case parameters.
void initialize() {
  running = false;
  draw_u = true;
  diffusion_only = false;
  constant_fk = true;
  min_u = 1;
  max_u = 1;
  min_v = 0;
  max_v = 0;
  println("(Re)Initialize");
  println("  Stop.");
  println("  Drawing u.");
  println("  Simulating reaction-diffusion.");
  println("  Using constant f, k:");
  println("    f = " + f);
  println("    k = " + k);
  for (int x = 0; x < cells_per_side; x++) {
    for (int y = 0; y < cells_per_side; y++) {
      u_curr[x][y] = 1;
      v_curr[x][y] = 0;
    }
  }
  
  // Vary u and v in a fixed rectangular region.
  for (int x = fixed_x; x < fixed_x + 10; x++) {
    for (int y = fixed_y; y < fixed_y + 10; y++) {
      float rand_u = .5 + random(-.05, .05);
      float rand_v = .25 + random(-.05, .05);
      u_curr[x][y] = rand_u;
      v_curr[x][y] = rand_v;
      
      // Update minmax.
      if (rand_u < min_u) {
        min_u = rand_u;
      } else if (rand_u > max_u) {
        max_u = rand_u;
      }
      if (rand_v < min_v) {
        min_v = rand_v;
      } else if (rand_v > max_v) {
        max_v= rand_v;
      }
    }
  }
  
  set_fk();
}

void set_fk() {
  for (int x = 0; x < cells_per_side; x++) {
    for (int y = 0; y < cells_per_side; y++) {
      if (constant_fk) {
        f_curr[x][y] = f;
        k_curr[x][y] = k;
      }
      else {
        float delta_f = (0 - .08)/(cells_per_side - 1);   // vary f from .08(top) -> 0(bottom)
        float delta_k = (.07 - .03)/(cells_per_side - 1); // vary k from .03(left) -> .07(right)
        f_curr[x][y] = .08 + (y * delta_f);         // decreasing top->bottom
        k_curr[x][y] = .03 + (x * delta_k);         // increasing left->right
      }
    }
  }
}

// Update u and v.
void reaction_diffusion() {
  
}

float[][] deep_copy(float[][] original) {
  float[][] copy = new float[original.length][original[0].length];
  for (int x = 0; x < original.length; x++) {
    for (int y = 0; y < original[0].length; y++) {
      copy[x][y] = original[x][y];
    }
  }
  return copy;
}

// Get the toroidal 8-neighbors of cell x, y.
int[][] get_toroidal_neighbors(int x, int y) {
  int[][] neighbors = new int[8][2];
  int n = 0;
  for (int xn = -1; xn <= 1; xn++) {
    for (int yn = -1; yn <= 1; yn++, n++) {
      if (abs(xn) == abs(yn)) {
        // if xn == yn == 0, it is the original cell.
        // if both xn and yn are either 1 or -1, then it is a diagonal neighbor.
        continue;
      }
      // Toroidal wrap of cell neighbors.
      int xx = (x + xn + cells_per_side) % cells_per_side;
      int yy = (y + yn + cells_per_side) % cells_per_side;
      
      neighbors[n] = new int[]{xx, yy};
    }
  }
  
  return neighbors;
}

void draw() {
  show_current_state();
  show_mouse_position();
  
  if (running) {
    reaction_diffusion(); 
  }
}

// Show u/v concentrations.
void show_current_state() {
  for (int x = 0; x < cells_per_side; x++) {
    for (int y = 0; y < cells_per_side; y++) {
      int gridX = cell_to_grid(x);
      int gridY = cell_to_grid(y);
      float min = min_u;
      float max = max_u;
      float curr = u_curr[x][y];
      if (!draw_u) {
        // Showing v concentration instead of u.
        min = min_v;
        max = max_v;
        curr = v_curr[x][y];
      }
      // Gradient min->max is black->white.
      float gradient = 255;
      if (min != max) gradient = ((curr - min)/(max - min)) * 255; // Avoid %0.
      fill(gradient);
      square(gridX, gridY, cell_size);
    }
  }
}

void show_mouse_position() {
  // Show location of mouse cell.
  int cellX = grid_to_cell(mouseX);
  int cellY = grid_to_cell(mouseY);
  int gridX = cell_to_grid(cellX);
  int gridY = cell_to_grid(cellY);
  stroke(red);
  noFill();
  square(gridX, gridY, cell_size);
  noStroke();
}

// Translate cell coordinates to grid coordinates.
int cell_to_grid(int cell) {
  return cell * cell_size;
}

// Translate grid coordinates to cell coordinates.
int grid_to_cell(int grid) {
  // Integer division results in the floor of the floating point value.
  return grid / cell_size; 
}

// Click on cells to turn them live.
void mousePressed() {
  int x = grid_to_cell(mouseX);
  int y = grid_to_cell(mouseY);
  println(x + ", " + y + ":");
  println("  u: " + u_curr[x][y]);
  println("  v: " + v_curr[x][y]);
  if (!constant_fk) {
    println("  f: " + f_curr[x][y]);
    println("  k: " + k_curr[x][y]);
  }
}

void keyPressed() {
  switch (key) {
    case('i'): 
      // Initialize the system with a fixed rectangular region that has specific u and v concentrations.
      initialize();
      break;
    case (' '):
      // Start or stop the simulation (toggle between these).
      println((!running? "Start." : "Stop") + "."); 
      running = !running;
      break;
    case ('u'):
      // At each timestep, draw values for u at each cell (default).
      println("Drawing u.");
      draw_u = true;
      break;
    case ('v'):
      // At each timestep, draw values for v at each cell (instead of u).
      println("Drawing v.");
      draw_u = false;
      break;
    case ('d'):
      // Toggle between performing diffusion alone or reaction-diffusion.
      println("Set simulation to " + (!diffusion_only? "diffusion only" : "reaction-diffusion") + "."); 
      diffusion_only = !diffusion_only;
      break;
    case ('p'):
      // Toggle between constant f, k for each cell and spatially-varying f, k.
      println("Set f, k to " + (!constant_fk? "constant: " : "varying.")); 
      if (!constant_fk) {
        println("  f = " + f);
        println("  k = " + k);
      }
      constant_fk = !constant_fk;
      set_fk();
      break;
    case ('1'):
      // Spots.
    case ('2'):
      // Stripes.
    case ('3'):
      // Spiral waves.
    case ('4'):
      // Custom pattern.
      
      // Change case.
      set_case(Character.getNumericValue(key) - 1);
      break;
    case ('t'):
      // Set boundary condition to toroidal.
      println("Set boundaries to toroidal.");
      boundary_condition = toroidal;
      break;
    case ('z'):
      // Set boundary condition to zero-derivative.
      println("Set boundaries to zero-derivative.");
      boundary_condition = zero_derivative;
      break;
    case ('f'):
      // Set boundary condition to fixed.
      println("Set boundaries to fixed.");
      boundary_condition = fixed;
      break;
    default:
      break;
  }
}
