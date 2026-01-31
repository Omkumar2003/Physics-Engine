import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

public class GravitySim3D {

    // --- SCREEN SETTINGS ---
    static int WIDTH = 400;  
    static int HEIGHT = 80;  
    static char[] SCREEN_BUFFER = new char[WIDTH * HEIGHT];
    static double[] Z_BUFFER = new double[WIDTH * HEIGHT]; 
    
    // --- PHYSICS SETTINGS ---
    static double G = 1.0;      
    static double DT = 0.02;     // Smaller time step for smoothness
    static int PHYSICS_SUBSTEPS = 4; // Run physics 4x per frame
    static double REBOUND = 0.8; 
    
    // --- RENDERING SETTINGS ---
    static char[] SHADE_CHARS = ".,-~:;=!*#$@".toCharArray();
    
    // Light coming from top-left-front
    static double LIGHT_X = 0.4;
    static double LIGHT_Y = -0.5;
    static double LIGHT_Z = 0.6; 

    static class Vec3 {
        double x, y, z;
        Vec3(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
        void add(double vx, double vy, double vz) { x += vx; y += vy; z += vz; }
    }

    static class Body {
        String name;
        Vec3 pos;
        Vec3 vel;
        double mass;
        double radius; 
        
        public Body(String name, double x, double y, double z, double vx, double vy, double vz, double mass, double radius) {
            this.name = name;
            this.pos = new Vec3(x, y, z);
            this.vel = new Vec3(vx, vy, vz);
            this.mass = mass;
            this.radius = radius;
        }

        public void update() {
            pos.x += vel.x * DT;
            pos.y += vel.y * DT;
            pos.z += vel.z * DT;
        }
    }

    static List<Body> bodies = new ArrayList<>();
    static Body sun;

    public static void main(String[] args) {
        // Normalize light vector
        double lMag = Math.sqrt(LIGHT_X*LIGHT_X + LIGHT_Y*LIGHT_Y + LIGHT_Z*LIGHT_Z);
        LIGHT_X /= lMag; LIGHT_Y /= lMag; LIGHT_Z /= lMag;

        // --- INIT SOLAR SYSTEM (FLAT DISK) ---
        
        // 1. The Sun
        sun = new Body("Sun", 0, 0, 0, 0, 0, 0, 10000, 10);
        bodies.add(sun);

        // 2. Mercury (Closest, Fast)
        // Positioned on X axis (28, 0, 0) -> Horizontal Orbit
        bodies.add(createOrbitingBody("Mercury", 28, 0, 0, 100, 2));

        // 3. Earth (Middle)
        // Positioned on Z axis (0, 0, 45) -> Horizontal Orbit
        bodies.add(createOrbitingBody("Earth", 0, 0, 45, 300, 4));

        // 4. Mars (Outer)
        // Positioned on negative X (-65, 0, 0) -> Horizontal Orbit
        bodies.add(createOrbitingBody("Mars", -65, 0, 0, 400, 5));
        
        // 5. Jupiter (Giant, Far out)
        // Positioned on negative Z (0, 0, -95) -> Horizontal Orbit
        bodies.add(createOrbitingBody("Jupiter", 0, 0, -95, 2000, 9));

        while (true) {
            for(int i=0; i<PHYSICS_SUBSTEPS; i++) physicsStep();
            renderStep();
            try { Thread.sleep(20); } catch (Exception e) {}
        }
    }

    // --- HELPER: AUTOMATIC HORIZONTAL ORBIT ---
    public static Body createOrbitingBody(String name, double x, double y, double z, double mass, double radius) {
        double dist = Math.sqrt(x*x + y*y + z*z);
        double vMag = Math.sqrt((G * sun.mass) / dist);

        // To make it horizontal (flat), we rotate the velocity 90 degrees around the Y axis.
        // If Pos is (X, 0, Z), Vel should be (Z, 0, -X) or (-Z, 0, X)
        
        double vx = z;
        double vy = 0; // Keep Y velocity 0 for flat orbits
        double vz = -x;

        // Normalize this direction vector
        double mag = Math.sqrt(vx*vx + vz*vz);
        vx /= mag;
        vz /= mag;

        // Apply speed
        vx *= vMag;
        vz *= vMag;

        return new Body(name, x, y, z, vx, vy, vz, mass, radius);
    }

    // --- PHYSICS ENGINE ---
    static void physicsStep() {
        for (int i = 0; i < bodies.size(); i++) {
            for (int j = i + 1; j < bodies.size(); j++) {
                Body b1 = bodies.get(i);
                Body b2 = bodies.get(j);

                double dx = b2.pos.x - b1.pos.x;
                double dy = b2.pos.y - b1.pos.y;
                double dz = b2.pos.z - b1.pos.z;
                
                double distSq = dx*dx + dy*dy + dz*dz;
                double dist = Math.sqrt(distSq);
                if (dist < 3) dist = 3; 

                double force = (G * b1.mass * b2.mass) / distSq;
                
                double fx = force * (dx / dist);
                double fy = force * (dy / dist);
                double fz = force * (dz / dist);

                b1.vel.add( fx / b1.mass * DT, fy / b1.mass * DT, fz / b1.mass * DT );
                b2.vel.add( -fx / b2.mass * DT, -fy / b2.mass * DT, -fz / b2.mass * DT );
            }
        }

        // Simple Collision
        for (int i = 0; i < bodies.size(); i++) {
            for (int j = i + 1; j < bodies.size(); j++) {
                Body b1 = bodies.get(i);
                Body b2 = bodies.get(j);
                double dx = b2.pos.x - b1.pos.x;
                double dy = b2.pos.y - b1.pos.y;
                double dz = b2.pos.z - b1.pos.z;
                double dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                double overlap = (b1.radius + b2.radius) - dist;
                
                if (overlap > 0) {
                   // Just separate them slightly to prevent stickiness
                   // (Full elastic bounce code removed for simplicity/performance in this view)
                   double nx = dx/dist; double ny = dy/dist; double nz = dz/dist;
                   b1.pos.add(-nx*overlap*0.5, -ny*overlap*0.5, -nz*overlap*0.5);
                   b2.pos.add(nx*overlap*0.5, ny*overlap*0.5, nz*overlap*0.5);
                }
            }
        }

        for (Body b : bodies) b.update();
    }

    // --- RENDER ENGINE ---
    static void renderStep() {
        Arrays.fill(SCREEN_BUFFER, ' ');
        Arrays.fill(Z_BUFFER, Double.MAX_VALUE);
        drawBorder();
        
        bodies.sort((a, b) -> Double.compare(b.pos.z, a.pos.z));

        for (Body b : bodies) renderSphere(b);

        System.out.print("\033[H");
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < HEIGHT; i++) {
            sb.append(SCREEN_BUFFER, i * WIDTH, WIDTH);
            sb.append('\n');
        }
        System.out.print(sb.toString());
    }

    static void renderSphere(Body b) {
        double focalLength = 80.0; 
        double zDepth = b.pos.z + 150; // Camera further back
        if (zDepth <= 1) return; 

        double scale = focalLength / zDepth;
        double screenCX = (b.pos.x * scale) * 2.0; 
        double screenCY = (b.pos.y * scale);
        double screenRad = b.radius * scale;
        
        int cx = (int)(WIDTH / 2 + screenCX);
        int cy = (int)(HEIGHT / 2 + screenCY);
        int rx = (int)Math.ceil(screenRad * 2.0);
        int ry = (int)Math.ceil(screenRad);

        int minX = Math.max(0, cx - rx);
        int maxX = Math.min(WIDTH - 1, cx + rx);
        int minY = Math.max(0, cy - ry);
        int maxY = Math.min(HEIGHT - 1, cy + ry);

        for (int y = minY; y <= maxY; y++) {
            for (int x = minX; x <= maxX; x++) {
                double ny = (y - cy) / (double)screenRad;
                double nx = (x - cx) / (double)(screenRad * 2.0);
                double distSq = nx*nx + ny*ny;
                
                if (distSq <= 1.0) {
                    double nz = Math.sqrt(1.0 - distSq);
                    double pixelDepth = zDepth - (nz * b.radius);
                    int idx = y * WIDTH + x;

                    if (idx >= 0 && idx < Z_BUFFER.length && pixelDepth < Z_BUFFER[idx]) {
                        Z_BUFFER[idx] = pixelDepth;
                        double dot = nx * LIGHT_X + ny * LIGHT_Y + nz * LIGHT_Z;
                        if (dot < 0) dot = 0;
                        int charIndex = (int)(dot * (SHADE_CHARS.length - 1));
                        SCREEN_BUFFER[idx] = SHADE_CHARS[Math.max(0, Math.min(charIndex, SHADE_CHARS.length-1))];
                    }
                }
            }
        }
    }

    static void drawBorder() {
        for(int x=0; x<WIDTH; x++) {
            SCREEN_BUFFER[x] = '-';
            SCREEN_BUFFER[(HEIGHT-1)*WIDTH + x] = '-';
        }
        for(int y=0; y<HEIGHT; y++) {
            SCREEN_BUFFER[y*WIDTH] = '|';
            SCREEN_BUFFER[y*WIDTH + (WIDTH-1)] = '|';
        }
    }
}
