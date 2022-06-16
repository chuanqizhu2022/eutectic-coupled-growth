#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <omp.h>
#include "CImg.h" //CImg ライブラリ（描画用）使用のためのヘッダ

using namespace std;
using namespace cimg_library;

#define N 3
#define NTH 8
#define NDX 128
#define NDY 128
#define NDZ 1
#define NDL 2560
#define PI 3.14159

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int mid = NDX / 2;
int rows = NDX / NTH;

int nstep = 10001;
int pstep = 500;

double dx = 1.0;
double dtime = 1.0;

double gamma0 = 0.1;
double astre = -0.05;
double mobi = 0.25;
double delta = 5.0 * dx;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double S1 = 0.03;
double S2 = 0.06;

double Dl = 0.1;
double Ds = 2.0e-4;

double gradT = 0.000;
double rateT = 0.000000;
double temp0 = -0.5;
double cl = 0.122;

double alpha_d = dtime * Dl / dx / dx;
double alpha_m = dtime / dx / dx * mobi * A0;

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];
double thij[N][N], vpij[N][N], etaij[N][N];
int anij[N][N];
double conlr[NDL], conlr2[NDL];

double
    ****phi,
    ****phi2,
    ***cont,
    ***cont2,
    ****conp,
    ***temp;
int
    ***phiNum,
    ****phiIdx;

int i, j, k, im, ip, jm, jp, km, kp;
int ni, phinum0;
int intpos, dist, hasS, allS, allL, fs1, fs2;
int curpos, prepos, frapass, intpass;
double invV;
double c0, c00, dc0, sum0, sumplane;
double sumc, suml, sums;

CImg<unsigned char> ch_fldxz(NDX, NDZ, 1, 3), ch_fldxy(NDX, NDY, 1, 3);
char outFileCh_xz[64], outFileCh_xy[64];

void datasave(int step);
double calC01e(double temp0), calC1e(double temp0), calC02e(double temp0), calC2e(double temp0);
double calDF10(double con0, double temp0, double dS), calDF20(double con0, double temp0, double dS);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "concenration field stablity number is: " << alpha_d << endl;
    cout << "phase field stability number is: " << alpha_m << endl;

    if ((alpha_d > 0.15) || (alpha_m > 0.15))
    {
        cout << "The computation is unstable, please change input parameters!" << endl;
        goto terminal;
    }

    // ---------------------------------  Initialization ------------------------------------
    phi = new double ***[N];
    phi2 = new double ***[N];
    conp = new double ***[N];
    for (ni = 0; ni <= nm; ni++)
    {
        phi[ni] = new double **[NDX];
        phi2[ni] = new double **[NDX];
        conp[ni] = new double **[NDX];
        for (i = 0; i <= ndmx; i++)
        {
            phi[ni][i] = new double *[NDY];
            phi2[ni][i] = new double *[NDY];
            conp[ni][i] = new double *[NDY];
            for (j = 0; j <= ndmy; j++)
            {
                phi[ni][i][j] = new double[NDZ];
                phi2[ni][i][j] = new double[NDZ];
                conp[ni][i][j] = new double[NDZ];
            }
        }
    }

    phiIdx = new int ***[N + 1];
    for (ni = 0; ni <= N; ni++)
    {
        phiIdx[ni] = new int **[NDX];
        for (i = 0; i <= ndmx; i++)
        {
            phiIdx[ni][i] = new int *[NDY];
            for (j = 0; j <= ndmy; j++)
            {
                phiIdx[ni][i][j] = new int[NDZ];
            }
        }
    }

    phiNum = new int **[NDX];
    for (i = 0; i <= ndmx; i++)
    {
        phiNum[i] = new int *[NDY];

        for (j = 0; j <= ndmy; j++)
        {
            phiNum[i][j] = new int[NDZ];
        }
    }

    cont = new double **[NDX];
    cont2 = new double **[NDX];
    temp = new double **[NDX];
    for (i = 0; i <= ndmx; i++)
    {
        cont[i] = new double *[NDY];
        cont2[i] = new double *[NDY];
        temp[i] = new double *[NDY];
        for (j = 0; j <= ndmy; j++)
        {
            cont[i][j] = new double[NDZ];
            cont2[i][j] = new double[NDZ];
            temp[i][j] = new double[NDZ];
        }
    }

    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            anij[i][j] = 0;
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
            }
        }
    }

    mij[1][2] = M0 * 0.1;
    mij[2][1] = M0 * 0.1;
    anij[1][0] = 1;
    anij[0][1] = 1;

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                temp[i][j][k] = temp0 + gradT * i * dx;
            }
        }
    }

    sum0 = 0.0;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                // if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 8 * NDX / 8)
                if (((i - NDX / 2) * (i - NDX / 2) + (k - NDZ / 2) * (k - NDZ / 2) > 400.0) && (j < NDY / 4))
                // if (i < NDX * 9.0 / 10.0 && j < NDY / 4)
                {
                    phi[1][i][j][k] = 1.0;
                    conp[1][i][j][k] = calC1e(temp[i][j][k]);
                    phi[2][i][j][k] = 0.0;
                    conp[2][i][j][k] = calC2e(temp[i][j][k]);
                    phi[0][i][j][k] = 0.0;
                    conp[0][i][j][k] = calC01e(temp[i][j][k]);
                }
                // else if (((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) >= (NDX * NDX / 2.0 / PI)) && (k < NDZ / 4))
                else if (((i - NDX / 2) * (i - NDX / 2) + (k - NDZ / 2) * (k - NDZ / 2) <= 400.0) && (j < NDY / 4))
                {
                    phi[1][i][j][k] = 0.0;
                    conp[1][i][j][k] = calC1e(temp[i][j][k]);
                    phi[2][i][j][k] = 1.0;
                    conp[2][i][j][k] = calC2e(temp[i][j][k]);
                    phi[0][i][j][k] = 0.0;
                    conp[0][i][j][k] = calC02e(temp[i][j][k]);
                }
                else
                {
                    phi[1][i][j][k] = 0.0;
                    conp[1][i][j][k] = calC1e(temp[i][j][k]);
                    phi[2][i][j][k] = 0.0;
                    conp[2][i][j][k] = calC2e(temp[i][j][k]);
                    phi[0][i][j][k] = 1.0;
                    conp[0][i][j][k] = cl;
                }
                cont[i][j][k] = conp[1][i][j][k] * phi[1][i][j][k] + conp[2][i][j][k] * phi[2][i][j][k] + conp[0][i][j][k] * phi[0][i][j][k];
                sum0 += cont[i][j][k];
            }
        }
    }
    c0 = sum0 / NDX / NDY / NDZ;
    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                if (i == ndmx)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndmx;
                }
                if (j == ndmy)
                {
                    jp = ndmy;
                }
                if (j == 0)
                {
                    jm = 0;
                }
                if (k == ndmz)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndmz;
                }

                phinum0 = 0;
                for (ni = 0; ni <= nm; ni++)
                {
                    if ((phi[ni][i][j][k] > 0.0) ||
                        ((phi[ni][i][j][k] == 0.0) && (phi[ni][ip][j][k] > 0.0) ||
                         (phi[ni][im][j][k] > 0.0) ||
                         (phi[ni][i][jp][k] > 0.0) ||
                         (phi[ni][i][jm][k] > 0.0) ||
                         (phi[ni][i][j][kp] > 0.0) ||
                         (phi[ni][i][j][km] > 0.0)))
                    {
                        phinum0++;
                        phiIdx[phinum0][i][j][k] = ni;
                    }
                }
                phiNum[i][j][k] = phinum0;
            }
        }
    }

#pragma omp parallel num_threads(NTH)
    {
        int istep, th_id;
        int start, end, offset;
        int startl, endL, offsetl;
        int ix, ixm, ixp, iy, iym, iyp, iz, izm, izp;
        int ii, jj, kk;
        int n1, n2, n3, phinum;
        double cddtt, sumcs, sumcl, con0izp;
        double dF, pddtt, psum, dsum;

        double phidx, phidy, phidz;
        double phidxx, phidyy, phidzz;
        double phidxy, phidxz, phidyz;
        double phiabs;

        double th, vp, eta;
        double epsilon0;

        double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
        double phidxp, phidyp, phidzp;
        double phidxpx, phidypx, phidzpx;
        double phidxpy, phidypy, phidzpy;
        double phidxpz, phidypz, phidzpz;

        double ep, epdx, epdy, epdz;

        double term0;
        double termx, termx0, termx1, termx0dx, termx1dx;
        double termy, termy0, termy1, termy0dy, termy1dy;
        double termz, termz0, termz1, termz0dz, termz1dz;

        double termiikk, termjjkk;

        istep = 0;
        th_id = omp_get_thread_num();

        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;

        startl = offsetl;

    start:;

        // ---------------------------------  Output the calculation data  ------------------------------------
        if ((istep % pstep == 0) && (th_id == 0))
        {
            intpass = frapass + intpos - NDZ / 4;
            curpos = intpass;
            invV = double(pstep) / double(curpos - prepos);
            prepos = curpos;
            datasave(istep);
            cout << "----------------------" << endl;
            cout << istep << " steps have done!" << endl;
            cout << "The interface position is " << intpos << endl;
            cout << "The average concnetration is " << c0 << endl;
            // ****** YZ *******
            cimg_forXY(ch_fldxz, x, z)
            {
                ch_fldxz(x, z, 0) = 255. * (cont[x][NDY / 2][z]); // red
                ch_fldxz(x, z, 1) = 255. * (cont[x][NDY / 2][z]); // green
                ch_fldxz(x, z, 2) = 255. * (cont[x][NDY / 2][z]); // blue
            }
            sprintf(outFileCh_xz, "figures/con_yz/2d%d.png", istep); // generate imagefile
            ch_fldxz.save_jpeg(outFileCh_xz);                        // save imagegilee

            // ****** XY *******
            cimg_forXY(ch_fldxy, x, y)
            {
                ch_fldxy(x, y, 0) = 255. * (cont[x][y][NDZ / 4]); // red
                ch_fldxy(x, y, 1) = 255. * (cont[x][y][NDZ / 4]); // green
                ch_fldxy(x, y, 2) = 255. * (cont[x][y][NDZ / 4]); // blue
            }
            sprintf(outFileCh_xy, "figures/con_xy/2d%d.png", istep); // generate imagefile
            ch_fldxy.save_jpeg(outFileCh_xy);                        // save imagegilee
        }

        // ---------------------------------  Evolution Equation of Phase fields ------------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
                    if (ix == ndmx)
                    {
                        ixp = 0;
                    }
                    if (ix == 0)
                    {
                        ixm = ndmx;
                    }
                    if (iy == ndmy)
                    {
                        iyp = ndmy;
                    }
                    if (iy == 0)
                    {
                        iym = 0;
                    }
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }

                    for (n1 = 1; n1 <= phiNum[ix][iy][iz]; n1++)
                    {
                        ii = phiIdx[n1][ix][iy][iz];
                        pddtt = 0.0;
                        for (n2 = 1; n2 <= phiNum[ix][iy][iz]; n2++)
                        {
                            jj = phiIdx[n2][ix][iy][iz];
                            dsum = 0.0;
                            for (n3 = 1; n3 <= phiNum[ix][iy][iz]; n3++)
                            {
                                kk = phiIdx[n3][ix][iy][iz];

                                phidx = (phi[kk][ixp][iy][iz] - phi[kk][ixm][iy][iz]) / 2.0;
                                phidy = (phi[kk][ix][iyp][iz] - phi[kk][ix][iym][iz]) / 2.0;
                                phidz = (phi[kk][ix][iy][izp] - phi[kk][ix][iy][izm]) / 2.0;

                                phidxx = (phi[kk][ixp][iy][iz] + phi[kk][ixm][iy][iz] - 2.0 * phi[kk][ix][iy][iz]);
                                phidyy = (phi[kk][ix][iyp][iz] + phi[kk][ix][iym][iz] - 2.0 * phi[kk][ix][iy][iz]);
                                phidzz = (phi[kk][ix][iy][izp] + phi[kk][ix][iy][izm] - 2.0 * phi[kk][ix][iy][iz]);

                                phidxy = (phi[kk][ixp][iyp][iz] + phi[kk][ixm][iym][iz] - phi[kk][ixm][iyp][iz] - phi[kk][ixp][iym][iz]) / 4.0;
                                phidxz = (phi[kk][ixp][iy][izp] + phi[kk][ixm][iy][izm] - phi[kk][ixm][iy][izp] - phi[kk][ixp][iy][izm]) / 4.0;
                                phidyz = (phi[kk][ix][iyp][izp] + phi[kk][ix][iym][izm] - phi[kk][ix][iym][izp] - phi[kk][ix][iyp][izm]) / 4.0;

                                phiabs = phidx * phidx + phidy * phidy + phidz * phidz;

                                if (anij[ii][kk] == 1 && phiabs != 0.0)
                                {
                                    epsilon0 = sqrt(aij[ii][kk]);

                                    th = thij[ii][kk];
                                    vp = vpij[ii][kk];
                                    eta = etaij[ii][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                    phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                    phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                    phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                    phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                    phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                    phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                    phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                    phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                    phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                    phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                    phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                    ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                    epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                    epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                    epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                    term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                    termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                    termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                    termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                    termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                    termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                    termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                    termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                    termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                    termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                    termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                    termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                    termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                    termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                    termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                    termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                    termiikk = term0 + termx + termy + termz;
                                }
                                else
                                {
                                    termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz) / (dx * dx);
                                }

                                if (anij[jj][kk] == 1 && phiabs != 0.0)
                                {
                                    epsilon0 = sqrt(aij[jj][kk]);

                                    th = thij[jj][kk];
                                    vp = vpij[jj][kk];
                                    eta = etaij[jj][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                    phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                    phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                    phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                    phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                    phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                    phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                    phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                    phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                    phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                    phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                    phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                    ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                    epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                    epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                    epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                    term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                    termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                    termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                    termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                    termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                    termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                    termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                    termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                    termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                    termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                    termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                    termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                    termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                    termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                    termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                    termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                    termjjkk = term0 + termx + termy + termz;
                                }
                                else
                                {
                                    termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz) / (dx * dx);
                                }

                                dsum += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][ix][iy][iz];
                            }
                            if (ii == 1 && jj == 0)
                            {
                                dF = calDF10(conp[0][ix][iy][iz], temp[ix][iy][iz], S1);
                            }
                            else if (ii == 0 && jj == 1)
                            {
                                dF = -calDF10(conp[0][ix][iy][iz], temp[ix][iy][iz], S1);
                            }
                            else if (ii == 2 && jj == 0)
                            {
                                dF = calDF20(conp[0][ix][iy][iz], temp[ix][iy][iz], S2);
                            }
                            else if (ii == 0 && jj == 2)
                            {
                                dF = -calDF20(conp[0][ix][iy][iz], temp[ix][iy][iz], S2);
                            }
                            else
                            {
                                dF = 0.0;
                            }
                            pddtt += -2.0 * mij[ii][jj] / double(phiNum[ix][iy][iz]) * (dsum - 8.0 / PI * dF * sqrt(phi[ii][ix][iy][iz] * phi[jj][ix][iy][iz]));
                        }
                        phi2[ii][ix][iy][iz] = phi[ii][ix][iy][iz] + pddtt * dtime;
                        if (phi2[ii][ix][iy][iz] >= 1.0)
                        {
                            phi2[ii][ix][iy][iz] = 1.0;
                        }
                        if (phi2[ii][ix][iy][iz] <= 0.0)
                        {
                            phi2[ii][ix][iy][iz] = 0.0;
                        }
                    }
                } // ix
            }     // iy
        }         // iz

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    for (kk = 0; kk <= nm; kk++)
                    {
                        phi[kk][ix][iy][iz] = phi2[kk][ix][iy][iz];
                    }
                }
            }
        }

        //
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    psum = 0.0;
                    for (kk = 0; kk <= nm; kk++)
                    {
                        psum += phi[kk][ix][iy][iz];
                    }
                    for (kk = 0; kk <= nm; kk++)
                    {
                        phi[kk][ix][iy][iz] = phi[kk][ix][iy][iz] / psum;
                        phi2[kk][ix][iy][iz] = phi[kk][ix][iy][iz];
                    }
                }
            }
        }
#pragma omp barrier
        // ---------------------------  Collect information of phase fields ----------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
                    if (ix == ndmx)
                    {
                        ixp = 0;
                    }
                    if (ix == 0)
                    {
                        ixm = ndmx;
                    }
                    if (iy == ndmy)
                    {
                        iyp = ndmy;
                    }
                    if (iy == 0)
                    {
                        iym = 0;
                    }
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }

                    phinum = 0;
                    for (ii = 0; ii <= nm; ii++)
                    {
                        if ((phi[ii][ix][iy][iz] > 0.0) ||
                            ((phi[ii][ix][iy][iz] == 0.0) &&
                                 (phi[ii][ixp][iy][iz] > 0.0) ||
                             (phi[ii][ixm][iy][iz] > 0.0) ||
                             (phi[ii][ix][iyp][iz] > 0.0) ||
                             (phi[ii][ix][iym][iz] > 0.0) ||
                             (phi[ii][ix][iy][izp] > 0.0) ||
                             (phi[ii][ix][iy][izm] > 0.0)))
                        {
                            phinum++;
                            phiIdx[phinum][ix][iy][iz] = ii;
                        }
                    }
                    phiNum[ix][iy][iz] = phinum;
                }
            }
        }
#pragma omp barrier
        // --------------------- Calculate concentration  --------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
                    if (ix == ndmx)
                    {
                        ixp = 0;
                    }
                    if (ix == 0)
                    {
                        ixm = ndmx;
                    }
                    if (iy == ndmy)
                    {
                        iyp = ndmy;
                    }
                    if (iy == 0)
                    {
                        iym = 0;
                    }
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }
                    if (phi[0][ix][iy][iz] == 0.0)
                    {
                        if (phi[1][ix][iy][iz] == 1.0)
                        {
                            conp[1][ix][iy][iz] = cont[ix][iy][iz];
                        }
                        else
                        {
                            conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                        }
                        if (phi[2][ix][iy][iz] == 1.0)
                        {
                            conp[2][ix][iy][iz] = cont[ix][iy][iz];
                        }
                        else
                        {
                            conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                        }
                        // Correct abnormal calculation at solid edge
                        if ((phi[0][ixp][iy][iz] > 0.0) || (phi[0][ixm][iy][iz] > 0.0) || (phi[0][ix][iyp][iz] > 0.0) || (phi[0][ix][iym][iz] > 0.0) || (phi[0][ix][iy][izp] > 0.0) || (phi[0][ix][iy][izm] > 0.0))
                        {
                            conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                            conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                        }
                        conp[0][ix][iy][iz] = calC01e(temp[ix][iy][iz]) * phi[1][ix][iy][iz] + calC02e(temp[ix][iy][iz]) * phi[2][ix][iy][iz];
                    }
                    else if (phi[0][ix][iy][iz] > 0.0 && phi[0][ix][iy][iz] < 1.0)
                    {
                        conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                        conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                        conp[0][ix][iy][iz] = (cont[ix][iy][iz] - conp[1][ix][iy][iz] * phi[1][ix][iy][iz] - conp[2][ix][iy][iz] * phi[2][ix][iy][iz]) / phi[0][ix][iy][iz];
                        // Correct abnormal calculation at liquid edge
                        if (phi[0][ix][iy][iz] < 0.05)
                        {
                            conp[0][ix][iy][iz] = (calC01e(temp[ix][iy][iz]) * phi[1][ix][iy][iz] + calC02e(temp[ix][iy][iz]) * phi[2][ix][iy][iz]) / (phi[1][ix][iy][iz] + phi[2][ix][iy][iz]);
                        }
                        if (conp[0][ix][iy][iz] > 1.0)
                        {
                            conp[0][ix][iy][iz] = 1.0;
                        }
                        if (conp[0][ix][iy][iz] < 0.0)
                        {
                            conp[0][ix][iy][iz] = 0.0;
                        }
                    }
                    else if (phi[0][ix][iy][iz] == 1.0)
                    {
                        conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                        conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                        conp[0][ix][iy][iz] = cont[ix][iy][iz];
                    }
                    cont[ix][iy][iz] = conp[1][ix][iy][iz] * phi[1][ix][iy][iz] + conp[2][ix][iy][iz] * phi[2][ix][iy][iz] + conp[0][ix][iy][iz] * phi[0][ix][iy][iz];
                    if (cont[ix][iy][iz] > 1.0)
                    {
                        cont[ix][iy][iz] = 1.0;
                    }
                    if (cont[ix][iy][iz] < 0.0)
                    {
                        cont[ix][iy][iz] = 0.0;
                    }
                }
            }
        }
        // --------------------- Correct concentration in liquid phase for mass conservation --------------------------
#pragma omp barrier
        // collect mass
        if (th_id == 0 && (istep > nstep / 5))
        {
            sumc = 0.0;
            suml = 0.0;
            for (ix = 0; ix <= ndmx; ix++)
            {
                for (iy = 0; iy <= ndmy; iy++)
                {
                    for (iz = 0; iz <= ndmz; iz++)
                    {
                        suml += phi[0][ix][iy][iz];
                        sumc += cont[ix][iy][iz];
                    }
                }
            }
            c00 = sumc / NDX / NDY / NDZ;
            dc0 = (sumc - c0 * NDX * NDY * NDZ) / suml;
            if (dc0 >= 0.001)
            {
                for (ix = 0; ix <= ndmx; ix++)
                {
                    for (iy = 0; iy <= ndmy; iy++)
                    {
                        for (iz = 0; iz <= ndmz; iz++)
                        {
                            conp[0][ix][iy][iz] = conp[0][ix][iy][iz] - dc0 * phi[0][ix][iy][iz];
                            cont[ix][iy][iz] = conp[0][ix][iy][iz] * phi[0][ix][iy][iz] + conp[1][ix][iy][iz] * phi[1][ix][iy][iz] + conp[2][ix][iy][iz] * phi[2][ix][iy][iz];
                            if (cont[ix][iy][iz] > 1.0)
                            {
                                cont[ix][iy][iz] = 1.0;
                            }
                            if (cont[ix][iy][iz] < 0.0)
                            {
                                cont[ix][iy][iz] = 0.0;
                            }
                        }
                    }
                }
            }
        }
#pragma omp barrier
        // ---------------------------------  Evolution Equation of Concentration field ------------------------------------
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    ixp = ix + 1;
                    ixm = ix - 1;
                    iyp = iy + 1;
                    iym = iy - 1;
                    izp = iz + 1;
                    izm = iz - 1;
                    if (ix == ndmx)
                    {
                        ixp = 0;
                    }
                    if (ix == 0)
                    {
                        ixm = ndmx;
                    }
                    if (iy == ndmy)
                    {
                        iyp = ndmy;
                    }
                    if (iy == 0)
                    {
                        iym = 0;
                    }
                    if (iz == ndmz)
                    {
                        izp = ndmz;
                    }
                    if (iz == 0)
                    {
                        izm = 0;
                    }
                    //拡散方程式内における微分計算
                    for (ii = 1; ii < nm; ii++)
                    {
                        sumcs = 0.25 * ((phi[ii][ixp][iy][iz] - phi[ii][ixm][iy][iz]) * (conp[ii][ixp][iy][iz] - conp[ii][ixm][iy][iz]) + (phi[ii][ix][iyp][iz] - phi[ii][ix][iym][iz]) * (conp[ii][ix][iyp][iz] - conp[ii][ix][iym][iz]) + (phi[ii][ix][iy][izp] - phi[ii][ix][iy][izm]) * (conp[ii][ix][iy][izp] - conp[ii][ix][iy][izm])) / dx / dx +
                                phi[ii][ix][iy][iz] * (conp[ii][ixp][iy][iz] + conp[ii][ixm][iy][iz] + conp[ii][ix][iyp][iz] + conp[ii][ix][iym][iz] + conp[ii][ix][iy][izp] + conp[ii][ix][iy][izm] - 6.0 * conp[ii][ix][iy][iz]) / dx / dx;
                    }
                    sumcl = 0.25 * ((phi[0][ixp][iy][iz] - phi[0][ixm][iy][iz]) * (conp[0][ixp][iy][iz] - conp[0][ixm][iy][iz]) + (phi[0][ix][iyp][iz] - phi[0][ix][iym][iz]) * (conp[0][ix][iyp][iz] - conp[0][ix][iym][iz]) + (phi[0][ix][iy][izp] - phi[0][ix][iy][izm]) * (conp[0][ix][iy][izp] - conp[0][ix][iy][izm])) / dx / dx +
                            phi[0][ix][iy][iz] * (conp[0][ixp][iy][iz] + conp[0][ixm][iy][iz] + conp[0][ix][iyp][iz] + conp[0][ix][iym][iz] + conp[0][ix][iy][izp] + conp[0][ix][iy][izm] - 6.0 * conp[0][ix][iy][iz]) / dx / dx;
                    cddtt = Ds * sumcs + Dl * sumcl;
                    cont2[ix][iy][iz] = cont[ix][iy][iz] + cddtt * dtime;
                    // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001; //濃度場の時間発展(陽解法)
                }
            }
        }

        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    cont[ix][iy][iz] = cont2[ix][iy][iz];
                }
            }
        }

        // cooling down
        for (ix = start; ix <= end; ix++)
        {
            for (iy = 0; iy <= ndmy; iy++)
            {
                for (iz = 0; iz <= ndmz; iz++)
                {
                    temp[ix][iy][iz] -= rateT * dtime;
                }
            }
        }

        //----------------------------------------------  Moving frame  -----------------------------------------------
#pragma omp barrier
        if (th_id == 0)
        {
            allL = 1;
            // search interface front
            for (iy = ndmy; iy >= 0; iy--)
            {
                if (allL == 0)
                {
                    intpos = iy;
                    break;
                }
                for (ix = 0; ix <= ndmx; ix++)
                {
                    for (iz = 0; iz <= ndmz; iz++)
                    {
                        if (phi[0][ix][iy][iz] < 1.0)
                        {
                            allL = 0;
                            break;
                        }
                    }
                    if (allL == 0)
                    {
                        break;
                    }
                }
            }
            istep = istep + 1;
            // check the distance from the middle of the domain
            if (intpos > mid)
            {
                // dist = intpos - mid;
                // frapass += dist;
                // cout << "--" << endl;
                // cout << "    the distance away from middle is " << dist << endl;
                // cout << "--" << endl;
                // cout << "    the interface temperature is " << temp[NDX / 2][NDY / 2][mid] << endl;

                // write passed domain
                // FILE *streamc;
                // char bufferc[30];
                // sprintf(bufferc, "data/con/passed.vtk");
                // streamc = fopen(bufferc, "a");

                // for (iz = 0; iz < dist; iz++)
                // {
                //     for (ix = 0; ix <= ndmx; ix++)
                //     {
                //         for (iy = 0; iy <= ndmy; iy++)
                //         {
                //             fprintf(streamc, "%e\n", cont[ix][iy][iz]);
                //         }
                //     }
                // }
                // fclose(streamc);

                for (iz = 0; iz <= (ndmz - dist); iz++)
                {
                    for (ix = 0; ix <= ndmx; ix++)
                    {
                        for (iy = 0; iy <= ndmy; iy++)
                        {
                            // temp
                            temp[ix][iy][iz] = temp[ix][iy][iz + dist];
                            // cont
                            cont[ix][iy][iz] = cont[ix][iy][iz + dist];
                            // phi
                            phi[0][ix][iy][iz] = phi[0][ix][iy][iz + dist];
                            phi[1][ix][iy][iz] = phi[1][ix][iy][iz + dist];
                            phi[2][ix][iy][iz] = phi[2][ix][iy][iz + dist];
                            // conp
                            conp[0][ix][iy][iz] = conp[0][ix][iy][iz + dist];
                            conp[1][ix][iy][iz] = conp[1][ix][iy][iz + dist];
                            conp[2][ix][iy][iz] = conp[2][ix][iy][iz + dist];
                        }
                    }
                }
                for (iz = (ndmz - dist + 1); iz <= ndmz; iz++)
                {
                    for (ix = 0; ix <= ndmx; ix++)
                    {
                        for (iy = 0; iy <= ndmy; iy++)
                        {
                            // temp
                            temp[ix][iy][iz] = temp[ix][iy][ndmz - dist] + gradT * (iz - ndmz + dist) * dx;
                            // cont
                            // new liquid is flowing into the box
                            cont[ix][iy][iz] = cl;
                            // phi
                            phi[0][ix][iy][iz] = 1.0;
                            phi[1][ix][iy][iz] = 0.0;
                            phi[2][ix][iy][iz] = 0.0;
                            // conp
                            conp[0][ix][iy][iz] = cl;
                            conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                            conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                        }
                    }
                }
            }
        }
#pragma omp barrier
        if (istep < nstep)
        {
            goto start;
        }
    end:;
    }
terminal:;
    return 0;
}

void datasave(int step)
{
    int i, j, k;

    FILE *streamc0;
    char bufferc0[30];
    sprintf(bufferc0, "data/con/3d%d.vtk", step);
    streamc0 = fopen(bufferc0, "a");

    fprintf(streamc0, "# vtk DataFile Version 1.0\n");
    fprintf(streamc0, "phi_%d.vtk\n", step);
    fprintf(streamc0, "ASCII\n");
    fprintf(streamc0, "DATASET STRUCTURED_POINTS\n");
    fprintf(streamc0, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
    fprintf(streamc0, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(streamc0, "ASPECT_RATIO 1.0 1.0 1.0\n");
    fprintf(streamc0, "\n");
    fprintf(streamc0, "POINT_DATA %d\n", NDX * NDY * NDZ);
    fprintf(streamc0, "SCALARS scalars float\n");
    fprintf(streamc0, "LOOKUP_TABLE default\n");

    for (k = 0; k <= ndmz; k++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (i = 0; i <= ndmx; i++)
            {
                // fprintf(streamc0, "%e\n", phi[1][i][j][k]);
                fprintf(streamc0, "%e\n", cont[i][j][k]);
            }
        }
    }
    fclose(streamc0);

    FILE *streamp;
    char bufferp[30];
    sprintf(bufferp, "data/phi/3d%d.vtk", step);
    streamp = fopen(bufferp, "a");

    fprintf(streamp, "# vtk DataFile Version 1.0\n");
    fprintf(streamp, "phi_%d.vtk\n", step);
    fprintf(streamp, "ASCII\n");
    fprintf(streamp, "DATASET STRUCTURED_POINTS\n");
    fprintf(streamp, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
    fprintf(streamp, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(streamp, "ASPECT_RATIO 1.0 1.0 1.0\n");
    fprintf(streamp, "\n");
    fprintf(streamp, "POINT_DATA %d\n", NDX * NDY * NDZ);
    fprintf(streamp, "SCALARS scalars float\n");
    fprintf(streamp, "LOOKUP_TABLE default\n");

    for (k = 0; k <= ndmz; k++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (i = 0; i <= ndmx; i++)
            {
                // fprintf(streamc0, "%e\n", phi[1][i][j][k]);
                fprintf(streamp, "%e\n", phi[2][i][j][k]);
            }
        }
    }
    fclose(streamp);
}

double calC01e(double temp0)
{
    double Te = 0.0;
    double ce = 0.122;
    double ml1 = -10.0;
    double kap1 = 0.131;
    double c01e;
    c01e = ce + (temp0 - Te) / ml1;
    return c01e;
}

double calC1e(double temp0)
{
    double Te = 0.0;
    double ce = 0.122;
    double ml1 = -10.0;
    double kap1 = 0.131;
    double c1e;
    c1e = (ce + (temp0 - Te) / ml1) * kap1;
    return c1e;
}

double calC02e(double temp0)
{
    double Te = 0.0;
    double ce = 0.122;
    double ml2 = 14.0;
    double c02e;
    c02e = ce + (temp0 - Te) / ml2;
    return c02e;
}

double calC2e(double temp0)
{
    double c2e = 1.0;
    return c2e;
}

double calDF10(double con0, double temp0, double dS)
{
    double Te = 0.0;
    double ce = 0.122;
    double ml1 = -10.0;
    double DF = ((con0 - ce) * ml1 + Te - temp0) * dS;
    return DF;
}

double calDF20(double con0, double temp0, double dS)
{
    double Te = 0.0;
    double ce = 0.122;
    double ml2 = 14.0;
    double DF = ((con0 - ce) * ml2 + Te - temp0) * dS;
    return DF;
}
