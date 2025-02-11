#include <iostream>
#include <vector>
#include <fstream> // Include for file handling
#include <iomanip> // Include for fixed precision formatting

//       ------->Y axis
//       |
//       |
//       |
//       V
//       x axis

using namespace std;

const int N = 10; // The size of the matrix (excluding boundary)

const double length = 10;

double h = length / N;

const double dt = 1;
const double diff = 0.1;

const double source = 100;
double a = diff * dt / (h * h);

vector<vector<double>> u(N + 2, vector<double>(N + 2));
vector<vector<double>> u0(N, vector<double>(N));

vector<vector<double>> v(N + 2, vector<double>(N + 2));
vector<vector<double>> v0(N + 2, vector<double>(N + 2));

vector<vector<double>> dens(N + 2, vector<double>(N + 2));
vector<vector<double>> dens0(N + 2, vector<double>(N + 2));

vector<vector<double>> xCo(N + 2, vector<double>(N + 2));
vector<vector<double>> yCo(N + 2, vector<double>(N + 2));

void SWAP(vector<vector<double>> &dens, vector<vector<double>> &dens0)
{
    vector<vector<double>> temp(N + 2, vector<double>(N + 2));
    temp = dens;
    dens0 = temp;
}

void createCoordinates(vector<vector<double>> &xCo, vector<vector<double>> &yCo)
{
    for (int i = 0; i < N + 1; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            xCo[i][j] = i * h;
            yCo[i][j] = j * h;
        }
    }
}
void addSource(vector<vector<double>> &dens)
{
    for (int i = 3; i <= 6; i++)
    {
        for (int j = 3; j <= 6; j++)
        {
            dens0[i][j] = source * dt;
        }
    }
}

void diffuse(vector<vector<double>> &dens, vector<vector<double>> &dens0)
{

    for (int k = 0; k < 300; k++)
    {
        for (int i = 1; i <= N; i++)
        {
            for (int j = 1; j <= N; j++)
            {
                dens[i][j] = (dens0[i][j] + a * (dens[i - 1][j] + dens[i + 1][j] + dens[i][j - 1] + dens[i][j + 1])) / (1 + 4 * a);
            }
        }
    }
}

int main()
{
    createCoordinates(xCo, yCo);

    // for (int i = 0; i < N + 1; i++)
    // {
    //     for (int j = 0; j < N + 1; j++)
    //     {
    //         cout << yCo[i][j] << " ";
    //     }
    //     cout << "\n";
    // }

    for (int t = 0; t < 10; t = t + dt)
    {
        addSource(dens);

        for (int i = 0; i <= N + 1; i++)
        {
            for (int j = 0; j <= N + 1; j++)
            {
                cout << dens0[i][j] << " ";
            }
            cout << "\n";
        }

        diffuse(dens, dens0);
        SWAP(dens, dens0);
        cout << "\n";
        cout << "\n";
    }
    cout << a;

    return 0;
}
