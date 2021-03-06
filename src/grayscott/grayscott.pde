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
// p/0 - Toggle between constant f, k for each cell (default) and spatially-varying f, k.
// 1 - Set parameters for spots (f = 0.035, k = 0.0625)
// 2 - Set parameters for stripes (f = 0.035, k = 0.06)
// 3 - Set parameters for spiral waves (f = 0.0118, k = 0.0475)
// 4-7 - Set parameters for custom patterns:
//   4 - concentric squares (f = .0271, k = .0557)
//   5 - flower (f = .0624, k = .0614)
//   6 - blinking spots (f = .0140, k = .0503)
//   7 - disappearing spots (f = .0108, k = .0531)
// (mouse click) - Print values for u and v at cell.  If in spatially-varying parameter mode, also print values for k and f at the cell.
//
// Custom Commands:
// t - Set boundary condition to periodic (toroidal, default).
// z - Set boundary condition to dirichlet (fixed-value).
// f - Set boundary condition to neumann (zero-derivative).
// >/. - Increase dt by .1 (up to 3; default 1)
// </, - Decrease dt by .1 (down to .1; default 1)


// Simulation Parameters.
float ru = 0.082;
float rv = 0.041;
float dt = 1;
float f = 0;
float k = 0;

String[] desc = new String[]{"spots", "stripes", "spiral waves", "custom - flower", 
                             "custom - honeycomb", "custom - blinking spots", "custom - disappearing spots"};
float[] f_case = new float[]{.035, .035, .0118, .0271, .0624, .0140, .0108};
float[] k_case = new float[]{.0625, .06, .0475, .0557, .0614, .0503, .0531};                             

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
int periodic = 0;
int dirichlet = 1;
int neumann = 2;
int boundary = periodic;    

// boundary constants for Dirichlet boundaries.
float boundary_constant_u = 1;
float boundary_constant_v = 0;

// Forward/implicit euler.
boolean implicit_euler = false;

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
  println("(Re)Initialize to case " + (curr_case + 1) + " (" + desc[curr_case] + ").");
  println("  Stop.");
  println("  Drawing u.");
  println("  Simulating reaction-diffusion.");
  println("  Using constant f, k:");
  println("    f = " + f);
  println("    k = " + k);
  print("  Boundary condition ");
  if (boundary == periodic) {
    println("periodic (toroidal).");
  }
  else if (boundary == neumann) {
    println("neumann (zero-derivative).");
  }
  else { // boundary == dirichlet
    println("dirichlet (fixed-value u = " + boundary_constant_u + ", v = " + boundary_constant_v + ").");
  }
  println("  dt = " + dt);
  for (int x = 0; x < cells_per_side; x++) {
    for (int y = 0; y < cells_per_side; y++) {
      u_curr[x][y] = 1 + random(-.05, .05);
      v_curr[x][y] = 0 + random(-.05, .05);
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
  u_next = deep_copy(u_curr);
  v_next = deep_copy(v_curr);
  // diffusion
  if (implicit_euler) {
    implicit_euler();
  }
  else {
    // forward Euler
    for (int x = 0; x < cells_per_side; x++) {
      for (int y = 0; y < cells_per_side; y++) {
        u_next[x][y] = u_curr[x][y] + dt * ru * Lu(x, y);
        v_next[x][y] = v_curr[x][y] + dt * rv * Lv(x, y);
      }
    }
  }
  if (!diffusion_only) {
    // reaction
    for (int x = 0; x < cells_per_side; x++) {
      for (int y = 0; y < cells_per_side; y++) {
        float uv2 = u_curr[x][y] * v_curr[x][y] * v_curr[x][y];
        float reaction_u = f_curr[x][y] * (1 - u_curr[x][y]) - uv2;
        float reaction_v = -(f_curr[x][y] + k_curr[x][y]) * v_curr[x][y] + uv2;
        u_next[x][y] = u_next[x][y] + dt *  reaction_u;
        v_next[x][y] = v_next[x][y] + dt * reaction_v;
      }
    }
  }
  
  u_curr = u_next;
  v_curr = v_next;
}

float Lu(int x, int y) {
  return L(x, y, u_curr, boundary_constant_u);
}

float Lv(int x, int y) {
  return L(x, y, v_curr, boundary_constant_v);
}

float L(int x, int y, float[][] vals, float boundary_constant) {
  int[][] neighbors = neighbors(x, y, boundary == periodic);
  float sum = 0;
  for (int i = 0; i < neighbors.length; i++) {
    int xx = neighbors[i][0];
    int yy = neighbors[i][1];
    sum += vals[xx][yy];
  }
  if (boundary == neumann) {
    // Ignore neighbors that don't exist
    return sum - neighbors.length * vals[x][y] ;
  }
  // If toroidal, neighbors.length is always 4.
  // If dirichlet, neighbors that don't exist are assigned a constant value.
  sum += boundary_constant * (4 - neighbors.length);
  return sum - 4 * vals[x][y] ;
}

void implicit_euler() {
  // TODO
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

// Get the toroidal 4-neighbors of cell x, y.
int[][] neighbors(int x, int y, boolean toroidal) {
  int[][] possible_neighbors = new int[4][2]; // Depending on the boundary condition, there may be < 4 neighbors.
  int n = 0;
  for (int xn = -1; xn <= 1; xn++) {
    for (int yn = -1; yn <= 1; yn++) {
      if (abs(xn) == abs(yn)) {
        // if xn == yn == 0, it is the original cell.
        // if both xn and yn are either 1 or -1, then it is a diagonal neighbor.
        continue;
      }
      // Toroidal wrap of cell neighbors.
      int xx = (x + xn + cells_per_side) % cells_per_side;
      int yy = (y + yn + cells_per_side) % cells_per_side;
      
      if (toroidal || ((abs(xx - x) == 1 || abs(xx - x) == 0) && (abs(yy - y) == 1 || abs(yy-y) == 0))) {
        // if xx - x != 1 or yy - y != 1, then the neighbor is toroidally wrapped.
        possible_neighbors[n] = new int[]{xx, yy};
        n++;
      }
    }
  }
  
  int[][] neighbors = new int[n][2];
  for (int i = 0; i < n; i++) {
    neighbors[i] = possible_neighbors[i];
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
    case ('0'):
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
      // Custom - flower
    case ('5'):
      // Custom - mesh.
    case ('6'):
      // Custom - blinking spots.
    case ('7'):
      // Custom - disappearing spots.
      
      // Change case.
      set_case(Character.getNumericValue(key) - 1);
      break;
    case ('t'):
      // Set boundary condition to periodic (toroidal).
      println("Set boundaries to periodic (toroidal).");
      boundary = periodic;
      break;
    case ('f'):
      // Set boundary condition to dirichlet (fixed-value).
      println("Set boundaries to dirichlet (fixed-value u = " + boundary_constant_u + ", v = " + boundary_constant_v + ").");
      boundary = dirichlet;
      break;
    case ('z'):
      // Set boundary condition to neumann (zero-derivative).
      println("Set boundaries to neumann (zero-derivative).");
      boundary = neumann;
      break;
    case ('e'):
      // Toggle between forward Euler and implicit Euler.
      println("Use " + (!implicit_euler? "implicit" : "forward") + " Euler."); 
      implicit_euler = !implicit_euler;
      break;
    case ('>'):
    case ('.'):
      // Increase dt by .1 (up to 3; default 1)
      dt = min(dt + .1, 3);
      println("dt = " + dt);
      break;
    case ('<'):
    case (','):
      // Decrease dt by .1 (down to .1; default 1)
      dt = max(dt - .1, .1);
      println("dt = " + dt);
      break;
    default:
      break;
  }
}
