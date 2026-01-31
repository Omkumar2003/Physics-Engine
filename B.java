import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

// ------------------- MAIN -------------------
public class B {
    public static void main(String[] args) throws InterruptedException {
        final int ROW = 120;
        final int COL = 100;

        Canvas canvas = new Canvas(ROW, COL);
        Physics physics = new Physics();

        List<Shape> shapes = new ArrayList<>();

        // Terrain
        Terrain terrain = new Terrain(ROW, COL);
        terrain.generate();
        shapes.add(terrain);

        // Balls
        Ball ball = new Ball(COL / 2.0, 0, 0.5, 0, 2, 0.8, physics.gravity);
        Ball ball1 = new Ball(COL / 4.0, 6, 1.5, 0, 3, 0.9, physics.gravity);

        shapes.add(ball);
        shapes.add(ball1);

        // Simulation loop
        for (int step = 0; step < 400; step++) {
            canvas.clear();

            // Update all shapes
            for (Shape s : shapes) {
                s.update(shapes, canvas, 0.1);
            }

            // Draw all shapes
            for (Shape s : shapes) {
                s.draw(canvas);
            }

            canvas.render();
            Thread.sleep(50);
        }
    }
}

// ------------------- CANVAS -------------------
class Canvas {
    int ROW, COL;
    char[][] pixels;

    public Canvas(int ROW, int COL) {
        this.ROW = ROW;
        this.COL = COL;
        pixels = new char[ROW][COL];
        clear();
    }

    public void clear() {
        for (int i = 0; i < ROW; i++) {
            for (int j = 0; j < COL; j++) {
                if (i == 0) pixels[i][j] = '-';
                else if (j == 0 || j == COL - 1) pixels[i][j] = '|';
                else pixels[i][j] = ' ';
            }
        }
    }

    public void render() {
        StringBuilder sb = new StringBuilder(ROW * (COL + 1));
        for (int i = 0; i < ROW; i++) {
            sb.append(pixels[i]);
            sb.append('\n');
        }
        System.out.print("\033[H");
        System.out.print(sb.toString());
        System.out.flush();
    }

    public void setPixel(int x, int y, char c) {
        if (x >= 0 && x < COL && y >= 0 && y < ROW) {
            pixels[y][x] = c;
        }
    }
}

// ------------------- PHYSICS -------------------
class Physics {
    double gravity;

    public Physics() {
        BigInteger earthMass = new BigInteger("5972200000000000000000000");
        BigInteger ballMass = new BigInteger("1");
        BigDecimal G = new BigDecimal("0.00000000006674");
        BigDecimal radiusOfEarth = new BigDecimal("6378137");
        BigDecimal radiusOfBall = new BigDecimal("1");
        BigDecimal r = radiusOfEarth.add(radiusOfBall);

        BigDecimal force = G.multiply(new BigDecimal(earthMass))
                .multiply(new BigDecimal(ballMass))
                .divide(r.pow(2), MathContext.DECIMAL128);
        BigDecimal acceleration = force.divide(new BigDecimal(ballMass), MathContext.DECIMAL128);
        gravity = acceleration.doubleValue();
    }
}

// ------------------- SHAPE -------------------
abstract class Shape {
    public abstract void update(List<Shape> shapes, Canvas canvas, double dt);
    public abstract void draw(Canvas canvas);
}

// ------------------- PHYSICAL SHAPE -------------------
abstract class PhysicalShape extends Shape {
    double x, y;
    double vX, vY;
    int radius;
    double damping;
    double gravity;
    double airResistance = 0.02;

    public PhysicalShape(double x, double y, double vX, double vY, int radius, double damping, double gravity) {
        this.x = x;
        this.y = y;
        this.vX = vX;
        this.vY = vY;
        this.radius = radius;
        this.damping = damping;
        this.gravity = gravity;
    }

    protected void applyPhysics(double dt) {
        vY += gravity * dt;
        vX -= vX * airResistance * dt;
        vY -= vY * airResistance * dt;
        x += vX * dt;
        y += vY * dt;
    }

    protected void checkWalls(Canvas canvas) {
        if (x + radius >= canvas.COL - 1) {
            x = canvas.COL - 1 - radius;
            vX = -vX * damping;
        }
        if (x - radius <= 0) {
            x = radius;
            vX = -vX * damping;
        }
        if (y - radius <= 0) {
            y = radius;
            vY = -vY * damping;
        }
    }

    protected void checkTerrainCollision(Terrain terrain) {
        int groundY = terrain.getHeight((int) x);
        if (y + radius >= groundY) {
            y = groundY - radius;
            vY = -vY * damping;
            vX += terrain.getSlope((int) x) * 0.05;
        }
    }
}

// ------------------- BALL -------------------
class Ball extends PhysicalShape {

    public Ball(double x, double y, double vX, double vY, int radius, double damping, double gravity) {
        super(x, y, vX, vY, radius, damping, gravity);
    }

    @Override
    public void update(List<Shape> shapes, Canvas canvas, double dt) {
        int subSteps = 5;  // higher = better collision resolution
        double dtSub = dt / subSteps;
        for (int i = 0; i < subSteps; i++) {
            applyPhysics(dtSub);

            for (Shape s : shapes) {
                if (s instanceof Terrain terrain) {
                    checkTerrainCollision(terrain);
                }
            }

            checkBallCollision(shapes);
            checkWalls(canvas);
        }
    }

    @Override
    public void draw(Canvas canvas) {
        int top = (int) y - radius;
        int bottom = (int) y + radius;
        if (top < 0) top = 0;
        if (bottom >= canvas.ROW) bottom = canvas.ROW - 1;

        for (int i = top; i <= bottom; i++) {
            double aspect = 2.5;
            int dx = (int) (Math.sqrt(radius * radius - (i - y) * (i - y)) * aspect);
            int left = (int) x - dx;
            int right = (int) x + dx;
            for (int j = left; j <= right; j++) {
                canvas.setPixel(j, i, '鬱');
            }
        }
    }
private void checkBallCollision(List<Shape> shapes) {
    for (Shape s : shapes) {
        if (s instanceof Ball other && other != this) {
            double dx = other.x - this.x;
            double dy = other.y - this.y;
            double distance = Math.sqrt(dx * dx + dy * dy);
            double minDist = this.radius + other.radius;

            if (distance < minDist && distance > 0) {
                // Normalize vector between balls
                double nx = dx / distance;
                double ny = dy / distance;

                // Separate the balls to prevent sticking
                double overlap = 0.5 * (minDist - distance);
                this.x -= overlap * nx;
                this.y -= overlap * ny;
                other.x += overlap * nx;
                other.y += overlap * ny;

                // Compute relative velocity along the normal
                double dvx = this.vX - other.vX;
                double dvy = this.vY - other.vY;
                double relVel = dvx * nx + dvy * ny;

                // Apply collision impulse (even if relative velocity is zero)
                double restitution = Math.min(this.damping, other.damping); // elasticity
                double impulse = -(1 + restitution) * relVel / 2; // equal mass
                this.vX += impulse * nx;
                this.vY += impulse * ny;
                other.vX -= impulse * nx;
                other.vY -= impulse * ny;
            }
        }
    }
}

}

// ------------------- TERRAIN -------------------
class Terrain extends Shape {
    int ROW, COL;
    int[] heights;
    Random rand = new Random();

    public Terrain(int ROW, int COL) {
        this.ROW = ROW;
        this.COL = COL;
        heights = new int[COL];
    }

    public void generate() {
        for (int j = 0; j < COL; j++) {
            int sine = (int) (Math.sin(j * 0.2) * 3);
            int bump = rand.nextInt(3);
            heights[j] = ROW - 5 - sine - bump;
        }
    }

    @Override
    public void update(List<Shape> shapes, Canvas canvas, double dt) {
        // Terrain does not move
    }

    @Override
    public void draw(Canvas canvas) {
        for (int j = 0; j < COL; j++) {
            canvas.setPixel(j, heights[j], '_');
        }
    }

    public int getHeight(int x) {
        if (x < 0) return heights[0];
        if (x >= COL) return heights[COL - 1];
        return heights[x];
    }

    public int getSlope(int x) {
        if (x < 0) x = 0;
        if (x >= COL - 1) x = COL - 2;
        return heights[x + 1] - heights[x];
    }
}
