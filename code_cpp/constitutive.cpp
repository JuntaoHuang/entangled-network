#include "constitutive.h"

double f(double lambda)
{

    return (lambda - 1);

    double lambda_lim = 5.0;
    double lambda_extension_ratio = 0.7;
    double lambda_max = lambda_lim * lambda_extension_ratio;
    if (lambda <= lambda_max)
    {
        return cohen_func(lambda, lambda_lim) - cohen_func(1.0, lambda_lim);
    }

    else
    {
        double f_lambda_m = 0.0;
        double df_lambda_m = 0.0;
        double d2f_lambda_m = 0.0;
        double d3f_lambda_m = 0.0;

        if (std::abs(lambda_extension_ratio - 0.9) < 1e-4)
        {
            f_lambda_m = 9.757017543860;
            df_lambda_m = 20.255401662050;
            d2f_lambda_m = 79.988336492200;
            d3f_lambda_m = 480.003683212989;
        }
        else if (std::abs(lambda_extension_ratio - 0.8) < 1e-4)
        {
            f_lambda_m = 4.627777777778;
            df_lambda_m = 5.261728395062;
            d2f_lambda_m = 9.986282578875;
            d3f_lambda_m = 30.004572473708;
        }
        else if (std::abs(lambda_extension_ratio - 0.7) < 1e-4)
        {
            f_lambda_m = 2.828431372549;
            df_lambda_m = 2.491426374471;
            d2f_lambda_m = 2.946679633022;
            d3f_lambda_m = 5.931672983552;
        }        
        else if (std::abs(lambda_extension_ratio - 0.6) < 1e-4)
        {
            f_lambda_m = 1.858333333333;
            df_lambda_m = 1.528125000000;
            d2f_lambda_m = 1.230468750000;
            d3f_lambda_m = 1.882324218750;
        }
        else if (std::abs(lambda_extension_ratio - 0.5) < 1e-4)
        {
            f_lambda_m = 1.216666666667;
            df_lambda_m = 1.088888888889;
            d2f_lambda_m = 0.616296296296;
            d3f_lambda_m = 0.777481481481;
        }
        else
        {
            std::cerr << "Error: lambda_extension_ratio is not supported in f(lambda)" << std::endl;
        }

        return f_lambda_m + df_lambda_m * (lambda - lambda_max) + 0.5 * d2f_lambda_m * std::pow(lambda - lambda_max, 2) + 1.0 / 6.0 * d3f_lambda_m * std::pow(lambda - lambda_max, 3);
    }
}

double energy_f(double lambda)
{

    return 0.5 * std::pow(lambda - 1, 2.0);

    double lambda_lim = 5.0;
    double lambda_extension_ratio = 0.7;
    double lambda_max = lambda_lim * lambda_extension_ratio;
    if (lambda <= lambda_max)
    {
        return relative_energy_cohen_func(lambda, lambda_lim) - relative_energy_cohen_func(1.0, lambda_lim);
    }
    else
    {
        double f_lambda_m = 0.0;
        double df_lambda_m = 0.0;
        double d2f_lambda_m = 0.0;
        double d3f_lambda_m = 0.0;
        double energy_correction = 0.0;

        if (std::abs(lambda_extension_ratio - 0.9) < 1e-4)
        {
            f_lambda_m = 9.757017543860;
            df_lambda_m = 20.255401662050;
            d2f_lambda_m = 79.988336492200;
            d3f_lambda_m = 480.003683212989;
            energy_correction = 7.866212728174;
        }
        else if (std::abs(lambda_extension_ratio - 0.8) < 1e-4)
        {
            f_lambda_m = 4.627777777778;
            df_lambda_m = 5.261728395062;
            d2f_lambda_m = 9.986282578875;
            d3f_lambda_m = 30.004572473708;
            energy_correction = 4.554146265059;
        }
        else if (std::abs(lambda_extension_ratio - 0.7) < 1e-4)
        {
            f_lambda_m = 2.828431372549;
            df_lambda_m = 2.491426374471;
            d2f_lambda_m = 2.946679633022;
            d3f_lambda_m = 5.931672983552;
            energy_correction = 2.745946127051;
        }        
        else if (std::abs(lambda_extension_ratio - 0.6) < 1e-4)
        {
            f_lambda_m = 1.858333333333;
            df_lambda_m = 1.528125000000;
            d2f_lambda_m = 1.230468750000;
            d3f_lambda_m = 1.882324218750;
            energy_correction = 1.593992207207;
        }
        else if (std::abs(lambda_extension_ratio - 0.5) < 1e-4)
        {
            f_lambda_m = 1.216666666667;
            df_lambda_m = 1.088888888889;
            d2f_lambda_m = 0.616296296296;
            d3f_lambda_m = 0.777481481481;
            energy_correction = 0.834300389658;
        }
        else
        {
            std::cerr << "Error: lambda_extension_ratio is not supported in energy_f(lambda)" << std::endl;
        }

        return f_lambda_m * (lambda - lambda_max) + 0.5 * df_lambda_m * std::pow(lambda - lambda_max, 2) + 1.0 / 6.0 * d2f_lambda_m * std::pow(lambda - lambda_max, 3) + 1.0 / 24.0 * d3f_lambda_m * std::pow(lambda - lambda_max, 4) + energy_correction;
    }    
}

double cohen_func(double x, double lambda_lim)
{
    return x / lambda_lim * (3.0 - std::pow(x / lambda_lim, 2.0)) / (1.0 - std::pow(x / lambda_lim, 2.0));
}

double relative_energy_cohen_func(double x, double lambda_lim)
{
    return std::pow(x, 2.0) / (2 * lambda_lim) - lambda_lim * std::log(lambda_lim * lambda_lim - x * x) - cohen_func(1, lambda_lim) * x;
}