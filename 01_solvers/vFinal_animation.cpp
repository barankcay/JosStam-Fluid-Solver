#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting
#include <cmath>

//       ------->Y axis
//       |
//       |
//       |
//       V
//       x axis

using namespace std;

const int N = 50; // The size of the matrix (excluding boundary)

const double length = 10;

double h = length / N;

const double dt = 1;
const double diff = 0.1;

const double visc = 0;
vector<vector<double>> u(N + 2, vector<double>(N + 2));
vector<vector<double>> u0(N, vector<double>(N));

vector<vector<double>> v(N + 2, vector<double>(N + 2));
vector<vector<double>> v0(N + 2, vector<double>(N + 2));

vector<vector<double>> dens(N + 2, vector<double>(N + 2));
vector<vector<double>> dens0(N + 2, vector<double>(N + 2));

vector<vector<double>> x(N + 2, vector<double>(N + 2));
vector<vector<double>> y(N + 2, vector<double>(N + 2));

vector<vector<double>> divergent(N + 2, vector<double>(N + 2));
vector<vector<double>> pressure(N + 2, vector<double>(N + 2));

void SWAP(vector<vector<double>> &dens, vector<vector<double>> &dens0)
{
    vector<vector<double>> temp(N + 2, vector<double>(N + 2));
    temp = dens;
    dens0 = temp;
}

void setBoundary(int b, vector<vector<double>> &M)
{
    for (int i = 1; i < N + 1; i++)
    {
        M[0][i] = b == 1 ? -M[1][i] : M[1][i];
        M[N + 1][i] = b == 1 ? -M[N][i] : M[N][i];
        M[i][0] = b == 2 ? -M[i][1] : M[i][1];
        M[i][N + 1] = b == 2 ? -M[i][N] : M[i][N];
    }
    M[0][0] = 0.5 * (M[1][0] + M[0][1]);
    M[0][N + 1] = 0.5 * (M[1][N + 1] + M[0][N]);
    M[N + 1][0] = 0.5 * (M[N][0] + M[N + 1][1]);
    M[N + 1][N + 1] = 0.5 * (M[N][N + 1] + M[N + 1][N]);
}

void velocInitialize(int xStart, int xEnd, int yStart, int yEnd, vector<vector<double>> &u, vector<vector<double>> &v, double uVeloc, double vVeloc)
{
    for (int i = xStart; i < xEnd; i++)
    {
        for (int j = yStart; j < yEnd; j++)
        {
            u[i][j] = uVeloc;
            v[i][j] = vVeloc;
        }
    }
}

double linearInterpolation(double a, double b, double c)
{
    double r;
    r = (c / h) * (b - a) + a;
    return r;
}

void createCoordinates(vector<vector<double>> &x, vector<vector<double>> &y)
{
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            x[i][j] = i * h;
            y[i][j] = j * h;
        }
    }
}
void addSource(int xStart, int xEnd, int yStart, int yEnd, vector<vector<double>> &dens, double source)
{
    for (int i = xStart; i < xEnd; i++)
    {
        for (int j = yStart; j < yEnd; j++)
        {
            dens[i][j] = source * dt;
        }
    }
}

void diffuse(int b, vector<vector<double>> &dens, vector<vector<double>> &dens0, double diff)
{
    double a = diff * dt / (h * h);

    for (int k = 0; k < 300; k++)
    {
        for (int i = 1; i < N + 1; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                dens[i][j] = (dens0[i][j] + a * (dens[i - 1][j] + dens[i + 1][j] + dens[i][j - 1] + dens[i][j + 1])) / (1 + 4 * a);
            }
        }
        setBoundary(b, dens);
    }
}

void advect(int b, vector<vector<double>> &dens, vector<vector<double>> &dens0, vector<vector<double>> &x, vector<vector<double>> &y, vector<vector<double>> &u, vector<vector<double>> &v)
{
    double x0;
    double y0;

    double z1;
    double z2;

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 1; j < N + 1; j++)
        {
            x0 = x[i][j] - u[i][j] * dt;
            y0 = y[i][j] - v[i][j] * dt;
            // cout << x0 << "\n";
            // cout << y0;

            // if (x0 < 0 || y0 < 0 || x0 > (N + 1) * h || y0 > (N + 1) * h)
            // {
            //     dens[i][j] = 0;
            //     continue;
            // }
            if (x0 < x[0][0])
            {
                x0 = x[0][0] + h / 2;
            }
            if (x0 > x[N + 1][0])
            {
                x0 = x[N + 1][0] - h / 2;
            }
            if (y0 < y[0][0])
            {
                y0 = y[0][0] + h / 2;
            }
            if (y0 > y[0][N + 1])
            {
                y0 = y[0][N + 1] - h / 2;
            }

            // z2 = dens0[ceil(x0 / h)][floor(y0 / h)] + (1 / h) * (y0 - y[ceil(x0 / h)][floor(y0 / h)]) * (dens0[ceil(x0 / h)][ceil(y0 / h)] - dens0[ceil(x0 / h)][floor(y0 / h)]);
            // z1 = dens0[floor(x0 / h)][floor(y0 / h)] + (1 / h) * (y0 - y[floor(x0 / h)][floor(y0 / h)]) * (dens0[floor(x0 / h)][ceil(y0 / h)] - dens0[floor(x0 / h)][floor(y0 / h)]);
            z2 = linearInterpolation(dens0[ceil(x0 / h)][floor(y0 / h)], dens0[ceil(x0 / h)][ceil(y0 / h)], y0 - y[ceil(x0 / h)][floor(y0 / h)]);
            z1 = linearInterpolation(dens0[floor(x0 / h)][floor(y0 / h)], dens0[floor(x0 / h)][ceil(y0 / h)], y0 - y[floor(x0 / h)][floor(y0 / h)]);
            dens[i][j] = linearInterpolation(z1, z2, x0 - x[floor(x0 / h)][floor(y0 / h)]);
            // cout << dens[i][j] << "\n";
        }
    }
    setBoundary(b, dens);
}

void project(int b, vector<vector<double>> &divergent, vector<vector<double>> &u, vector<vector<double>> &v, vector<vector<double>> &pressure)
{
    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 1; j < N + 1; j++)
        {
            divergent[i][j] = (u[i + 1][j] - u[i - 1][j] + v[i][j + 1] - v[i][j - 1]) / 2;
        }
    }
    setBoundary(b, divergent);
    setBoundary(b, pressure);
    for (int k = 0; k < 300; k++)
    {
        for (int i = 1; i < N + 1; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                pressure[i][j] = (pressure[i - 1][j] + pressure[i + 1][j] + pressure[i][j - 1] + pressure[i][j + 1] - divergent[i][j]) / 4;
            }
        }
        setBoundary(b, pressure);
    }

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 1; j < N + 1; j++)
        {
            u[i][j] -= 0.5 * (pressure[i + 1][j] - pressure[i - 1][j]);
            v[i][j] -= 0.5 * (pressure[i][j + 1] - pressure[i][j - 1]);
        }
    }
    setBoundary(2, u);
    setBoundary(1, v);
}

void density_step(vector<vector<double>> &dens, vector<vector<double>> &dens0, vector<vector<double>> &u, vector<vector<double>> &v, vector<vector<double>> &x, vector<vector<double>> &y)
{
    // Four lines below are for adding source
    int xStart = (N / 2) - 3;
    int xEnd = (N / 2) + 3;
    int yStart = (N / 2) - 3;
    int yEnd = (N / 2) + 3;

    double source = 300;
    addSource(xStart, xEnd, yStart, yEnd, dens, source);
    SWAP(dens, dens0);
    diffuse(0, dens, dens0, diff);

    SWAP(dens, dens0);
    advect(0, dens, dens0, x, y, u, v);

    addSource((N / 2) - 3, (N / 2) + 3, (N / 2) - 3, (N / 2) + 3, dens, 300);
}

void velocity_step(vector<vector<double>> &u, vector<vector<double>> &v, vector<vector<double>> &u0, vector<vector<double>> &v0, vector<vector<double>> &divergent, vector<vector<double>> &pressure, vector<vector<double>> &x, vector<vector<double>> &y)
{
    // Six lines below are for initializing the velocity field
    int xStart = (N / 2) - 3;
    int xEnd = (N / 2) + 3;
    int yStart = 3;
    int yEnd = 4;
    double xVeloc = 0;
    double yVeloc = 3;

    velocInitialize(xStart, xEnd, yStart, yEnd, u, v, xVeloc, yVeloc);

    SWAP(u, u0);
    diffuse(2, u, u0, visc);

    SWAP(v, v0);
    diffuse(1, v, v0, visc);

    project(0, divergent, u, v, pressure);
    SWAP(u, u0);
    SWAP(v, v0);

    advect(2, u, u0, x, y, u, v);

    advect(1, v, v0, x, y, u, v);

    project(0, divergent, u, v, pressure);
}

// Function to save dens0 to a file as a matrix in CSV format for Excel
// Function to save dens0 to a file as a matrix in CSV format for Excel
void saveToFile(const vector<vector<double>> &dens, const string &filename)
{
    ofstream outFile(filename);
    if (outFile.is_open())
    {
        // Set fixed point notation and set precision for writing to the file
        outFile << fixed << setprecision(6); // Set the precision to 6 decimal places

        // Write the data row by row, each row being a line in the CSV
        for (int i = 0; i <= N + 1; i++) // Include boundary cells (0 to N+1)
        {
            for (int j = 0; j <= N + 1; j++) // Include boundary cells (0 to N+1)
            {
                outFile << dens[i][j]; // Write the value

                if (j < N + 1)      // Avoid adding a comma at the end of the row
                    outFile << ","; // Separate values with a comma
            }
            outFile << "\n"; // New line for each row
        }

        outFile.close();
        cout << "Data saved to " << filename << endl;
    }
    else
    {
        cerr << "Unable to open file: " << filename << endl;
    }
}

int main()
{
    createCoordinates(x, y);

    for (int t = 0; t < 300; t = t + dt)
    {
        saveToFile(dens, "dens_t" + to_string(t) + ".csv");
        saveToFile(v, "v_t" + to_string(t) + ".csv");
        velocity_step(u, v, u0, v0, divergent, pressure, x, y);
        density_step(dens, dens0, u, v, x, y);

        // for (int i = 0; i <= N + 1; i++) // Include boundary cells (0 to N+1)
        // {
        //     for (int j = 0; j <= N + 1; j++) // Include boundary cells (0 to N+1)
        //     {
        //         cout << u[i][j] << " ";
        //     }
        //     cout << "\n";
        // }

        // cout << "\n";
        // cout << "\n";

        // Save the current dens0 to a file after each time step
    }

    return 0;
}
