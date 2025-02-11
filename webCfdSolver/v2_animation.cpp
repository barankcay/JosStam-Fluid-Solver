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

const int N = 20; // The size of the matrix (excluding boundary)

const double length = 10;

double h = length / N;

const double dt = 1;
const double diff = 0.1;

double a = diff * dt / (h * h);

vector<vector<double>> u(N + 2, vector<double>(N + 2));
vector<vector<double>> u0(N, vector<double>(N));

vector<vector<double>> v(N + 2, vector<double>(N + 2));
vector<vector<double>> v0(N + 2, vector<double>(N + 2));

vector<vector<double>> dens(N + 2, vector<double>(N + 2));
vector<vector<double>> dens0(N + 2, vector<double>(N + 2));

vector<vector<double>> x(N + 2, vector<double>(N + 2));
vector<vector<double>> y(N + 2, vector<double>(N + 2));

void SWAP(vector<vector<double>> &dens, vector<vector<double>> &dens0)
{
    vector<vector<double>> temp(N + 2, vector<double>(N + 2));
    temp = dens;
    dens0 = temp;
}

void velocInitialize(vector<vector<double>> &u, vector<vector<double>> &v, double uVeloc, double vVeloc)
{
    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 1; j < N + 1; j++)
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
void addSource(int xStart, int xEnd, int yStart, int yEnd, vector<vector<double>> &dens0, double source)
{
    for (int i = xStart; i < xEnd; i++)
    {
        for (int j = yStart; j < yEnd; j++)
        {
            dens0[i][j] = source * dt;
        }
    }
}

void diffuse(vector<vector<double>> &dens, vector<vector<double>> &dens0)
{

    for (int k = 0; k < 300; k++)
    {
        for (int i = 1; i < N + 1; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                dens[i][j] = (dens0[i][j] + a * (dens[i - 1][j] + dens[i + 1][j] + dens[i][j - 1] + dens[i][j + 1])) / (1 + 4 * a);
            }
        }
    }
}

void advect(vector<vector<double>> &dens, vector<vector<double>> &dens0, vector<vector<double>> &x, vector<vector<double>> &y, vector<vector<double>> &u, vector<vector<double>> &v)
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

            if (x0 <= 0 || y0 <= 0 || x0 > (N + 1) * h || y0 > (N + 1) * h)
            {
                dens[i][j] = 0;
                continue;
            }

            // z2 = dens0[ceil(x0 / h)][floor(y0 / h)] + (1 / h) * (y0 - y[ceil(x0 / h)][floor(y0 / h)]) * (dens0[ceil(x0 / h)][ceil(y0 / h)] - dens0[ceil(x0 / h)][floor(y0 / h)]);
            // z1 = dens0[floor(x0 / h)][floor(y0 / h)] + (1 / h) * (y0 - y[floor(x0 / h)][floor(y0 / h)]) * (dens0[floor(x0 / h)][ceil(y0 / h)] - dens0[floor(x0 / h)][floor(y0 / h)]);
            z2 = linearInterpolation(dens0[ceil(x0 / h)][floor(y0 / h)], dens0[ceil(x0 / h)][ceil(y0 / h)], y0 - y[ceil(x0 / h)][floor(y0 / h)]);
            z1 = linearInterpolation(dens0[floor(x0 / h)][floor(y0 / h)], dens0[floor(x0 / h)][ceil(y0 / h)], y0 - y[floor(x0 / h)][floor(y0 / h)]);
            dens[i][j] = linearInterpolation(z1, z2, x0 - x[floor(x0 / h)][floor(y0 / h)]);
            // cout << dens[i][j] << "\n";
        }
    }
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
    velocInitialize(u, v, 0, 0.3);
    addSource(0, 19, 0, 19, dens0, 100);

    for (int t = 0; t < 100; t = t + dt)
    {

        saveToFile(dens0, "dens0_t" + to_string(t) + ".csv");
        // Print dens0 in the console for each timestep (optional)
        for (int i = 0; i <= N + 1; i++) // Include boundary cells (0 to N+1)
        {
            for (int j = 0; j <= N + 1; j++) // Include boundary cells (0 to N+1)
            {
                cout << dens0[i][j] << " ";
            }
            cout << "\n";
        }

        diffuse(dens, dens0);
        SWAP(dens, dens0);
        advect(dens, dens0, x, y, u, v);
        SWAP(dens, dens0);
        cout << "\n";
        cout << "\n";

        // Save the current dens0 to a file after each time step
    }
    cout << a;

    return 0;
}
