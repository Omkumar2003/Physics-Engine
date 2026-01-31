import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class FluidSPH {

    // --- CONFIGURATION ---
    // Console Dimensions (Adjust to fit your terminal)
    static final int WIDTH = 300;
    static final int HEIGHT = 100;

    // Simulation Constants (Tuned for ASCII "Physics")
    static final double GRAVITY = 0.15;
    static final double SMOOTHING_RADIUS = 5.0; // The "h" in the math
    static final double TARGET_DENSITY = 0.8;
    static final double PRESSURE_MULTIPLIER = 80.0; // The stiffness
    static final double NEAR_PRESSURE_MULTIPLIER = 80.0; 
    static final double VISCOSITY_STRENGTH = 0.05;
    static final double DAMPING = 0.5; // Energy loss on wall collision

    // Time step
    static final double DT = 0.5;

    public static void main(String[] args) throws InterruptedException {
        FluidSystem system = new FluidSystem();
        ConsoleRenderer renderer = new ConsoleRenderer(WIDTH, HEIGHT);

        // ---------------------------------------------------------
        // SCENARIO: DAM BREAK
        // Create a block of water particles on the left side
        // ---------------------------------------------------------
        int cols = 20;
        int rows = 30;
        double startX = 5;
        double startY = HEIGHT - 2;

        for (int y = 0; y < rows; y++) {
            for (int x = 0; x < cols; x++) {
                // Add slight jitter to positions to prevent perfect stacking glitches
                double jitter = Math.random() * 0.5;
                system.addParticle(startX + x * 1.5 + jitter, startY - y * 1.2);
            }
        }

        System.out.println("Simulation Started. Press Ctrl+C to stop.");
        Thread.sleep(1000);

        // --- MAIN LOOP ---
        while (true) {
            long start = System.currentTimeMillis();

            // 1. Update Physics
            system.update();

            // 2. Render
            renderer.clear();
            renderer.drawParticles(system.particles);
            renderer.render();

            // 3. FPS Cap (Approx 30 FPS)
            long elapsed = System.currentTimeMillis() - start;
            long wait = 33 - elapsed;
            if (wait > 0) Thread.sleep(wait);
        }
    }
}

// ---------------------------------------------------------
// DATA STRUCTURES
// ---------------------------------------------------------

class Particle {
    double x, y;
    double vx, vy; // Velocity
    double fx, fy; // Force
    double density;
    double nearDensity;
    double pressure;
    double nearPressure;

    public Particle(double x, double y) {
        this.x = x;
        this.y = y;
    }
}

// ---------------------------------------------------------
// PHYSICS ENGINE (Smoothed Particle Hydrodynamics)
// ---------------------------------------------------------

class FluidSystem {
    List<Particle> particles = new ArrayList<>();

    public void addParticle(double x, double y) {
        particles.add(new Particle(x, y));
    }

    public void update() {
        // Step 1: Calculate Density & Pressure
        // (Corresponds to timestamp 10:47 in the video)
        computeDensityPressure();

        // Step 2: Calculate Internal Forces (Pressure + Viscosity)
        // (Corresponds to timestamp 17:00 and 34:00)
        computeForces();

        // Step 3: Integrate (Move particles) & Enforce Boundaries
        integrate();
    }

    private void computeDensityPressure() {
        for (Particle p : particles) {
            p.density = 0;
            p.nearDensity = 0;

            for (Particle other : particles) {
                if (p == other) continue; // Skip self

                double dist = Math.hypot(p.x - other.x, p.y - other.y);
                
                if (dist < FluidSPH.SMOOTHING_RADIUS) {
                    // Standard Density
                    p.density += smoothingKernel(dist, FluidSPH.SMOOTHING_RADIUS);
                    
                    // Near Density (Simulating "surface tension" / preventing clustering)
                    // Timestamp 37:54 in video
                    p.nearDensity += nearSmoothingKernel(dist, FluidSPH.SMOOTHING_RADIUS);
                }
            }

            // Convert Density to Pressure
            // P = k * (density - target_density)
            p.pressure = FluidSPH.PRESSURE_MULTIPLIER * (p.density - FluidSPH.TARGET_DENSITY);
            p.nearPressure = FluidSPH.NEAR_PRESSURE_MULTIPLIER * p.nearDensity;
        }
    }

    private void computeForces() {
        for (Particle p : particles) {
            p.fx = 0;
            p.fy = 0;
            
            // Gravity
            p.fy += FluidSPH.GRAVITY;

            for (Particle other : particles) {
                if (p == other) continue;

                double dx = other.x - p.x;
                double dy = other.y - p.y;
                double dist = Math.sqrt(dx * dx + dy * dy);

                if (dist < FluidSPH.SMOOTHING_RADIUS && dist > 0.0001) {
                    // Normalized direction vector
                    double dirX = dx / dist;
                    double dirY = dy / dist;

                    // --- PRESSURE FORCE ---
                    // Calculate shared pressure (Newton's 3rd Law optimization)
                    // Timestamp 21:40
                    double sharedPressure = (p.pressure + other.pressure) / 2.0;
                    double sharedNearPressure = (p.nearPressure + other.nearPressure) / 2.0;

                    // Apply Force based on Slope of kernels
                    double slope = smoothingKernelDerivative(dist, FluidSPH.SMOOTHING_RADIUS);
                    double nearSlope = nearSmoothingKernelDerivative(dist, FluidSPH.SMOOTHING_RADIUS);

                    double forceMagnitude = (sharedPressure * slope) + (sharedNearPressure * nearSlope);

                    // Apply to total force (F = ma, assume m=1)
                    // We divide by density to convert pressure to acceleration
                    if (other.density > 0.001) {
                        double acc = forceMagnitude / other.density;
                        p.fx += dirX * acc;
                        p.fy += dirY * acc;
                    }

                    // --- VISCOSITY FORCE ---
                    // Timestamp 34:42
                    // Smoothes relative velocities
                    double viscWeight = smoothingKernel(dist, FluidSPH.SMOOTHING_RADIUS);
                    double vxDiff = other.vx - p.vx;
                    double vyDiff = other.vy - p.vy;

                    p.fx += vxDiff * FluidSPH.VISCOSITY_STRENGTH * viscWeight;
                    p.fy += vyDiff * FluidSPH.VISCOSITY_STRENGTH * viscWeight;
                }
            }
        }
    }

    private void integrate() {
        for (Particle p : particles) {
            // Update Velocity
            p.vx += p.fx * FluidSPH.DT;
            p.vy += p.fy * FluidSPH.DT;

            // Update Position
            p.x += p.vx * FluidSPH.DT;
            p.y += p.vy * FluidSPH.DT;

            // --- BOUNDARY CHECKS ---
            // Floor
            if (p.y >= FluidSPH.HEIGHT - 1) {
                p.y = FluidSPH.HEIGHT - 1.01;
                p.vy *= -FluidSPH.DAMPING;
                p.vx *= 0.9; // Friction with floor
            }
            // Ceiling
            else if (p.y <= 1) {
                p.y = 1.01;
                p.vy *= -FluidSPH.DAMPING;
            }

            // Walls
            if (p.x <= 1) {
                p.x = 1.01;
                p.vx *= -FluidSPH.DAMPING;
            } else if (p.x >= FluidSPH.WIDTH - 1) {
                p.x = FluidSPH.WIDTH - 1.01;
                p.vx *= -FluidSPH.DAMPING;
            }
        }
    }

    // --- SPH KERNELS ---
    // These define the "shape" of the influence a particle has on neighbors

    // Poly6 Kernel (Standard Smoothing)
    private double smoothingKernel(double dist, double radius) {
        if (dist >= radius) return 0;
        double volume = Math.PI * Math.pow(radius, 4) / 6;
        return Math.pow(radius * radius - dist * dist, 2) / volume;
    }

    // Spiky Kernel Derivative (For Pressure Gradient)
    // Returns negative value (pushing away)
    private double smoothingKernelDerivative(double dist, double radius) {
        if (dist >= radius) return 0;
        double scale = 12 / (Math.pow(radius, 4) * Math.PI);
        return -(dist - radius) * scale;
    }

    // Near-Density Kernel (Spikier to prevent overlap)
    private double nearSmoothingKernel(double dist, double radius) {
        if (dist >= radius) return 0;
        double v = 1 - (dist / radius);
        return v * v * v;
    }

    private double nearSmoothingKernelDerivative(double dist, double radius) {
        if (dist >= radius) return 0;
        double v = 1 - (dist / radius);
        return -3 * v * v / radius;
    }
}

// ---------------------------------------------------------
// RENDERER (ANSI ART)
// ---------------------------------------------------------

class ConsoleRenderer {
    int w, h;
    char[][] buffer;
    String[][] colorBuffer;

    // ANSI Colors
    static final String ANSI_RESET = "\u001B[0m";
    static final String ANSI_BLUE = "\u001B[34m";
    static final String ANSI_CYAN = "\u001B[36m";
    static final String ANSI_WHITE = "\u001B[37m";
    static final String ANSI_RED = "\u001B[31m";
    static final String ANSI_GRAY = "\u001B[90m";

    public ConsoleRenderer(int w, int h) {
        this.w = w;
        this.h = h;
        buffer = new char[h][w];
        colorBuffer = new String[h][w];
    }

    public void clear() {
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (x == 0 || x == w - 1 || y == 0 || y == h - 1) {
                    buffer[y][x] = '#';
                    colorBuffer[y][x] = ANSI_GRAY;
                } else {
                    buffer[y][x] = ' ';
                    colorBuffer[y][x] = ANSI_RESET;
                }
            }
        }
    }

    public void drawParticles(List<Particle> particles) {
        for (Particle p : particles) {
            int ix = (int) p.x;
            int iy = (int) p.y;

            if (ix > 0 && ix < w - 1 && iy > 0 && iy < h - 1) {
                // Visualisation: Change color based on pressure
                // High pressure (compressed) = White/Red
                // Low pressure (relaxed) = Blue
                
                String color;
                char symbol;

                if (p.pressure > 60) {
                    color = ANSI_RED;
                    symbol = '@'; // Highly compressed
                } else if (p.pressure > 30) {
                    color = ANSI_WHITE;
                    symbol = '0';
                } else {
                    color = ANSI_BLUE;
                    symbol = 'o';
                }

                // Simple depth buffer check (if multiple particles in one cell, show the high pressure one)
                if (buffer[iy][ix] == ' ' || buffer[iy][ix] == 'o') {
                    buffer[iy][ix] = symbol;
                    colorBuffer[iy][ix] = color;
                }
            }
        }
    }

    public void render() {
        StringBuilder sb = new StringBuilder(w * h + 20);
        // Reset cursor to top-left
        sb.append("\033[H");

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                sb.append(colorBuffer[y][x]).append(buffer[y][x]);
            }
            sb.append("\n");
        }
        sb.append(ANSI_RESET);
        System.out.print(sb.toString());
        System.out.flush();
    }
}
