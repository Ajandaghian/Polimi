#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.diagnostic import acorr_ljungbox
import statsmodels.api as sm


def IMR(original_df, col_name):
    """Implements the Individual Moving Range (IMR) chart.
    Parameters
    ----------
    original_df : pandas.DataFrame
        A DataFrame containing the data to be plotted.
    col_name : str
        The name of the column to be used for the IMR control chart.
    Returns
    -------
    chart : matplotlib.axes._subplots.AxesSubplot
        The IMR chart.

    df_IMR : pandas.DataFrame with the following additional columns
        - MR: The moving range
        - I_UCL: The upper control limit for the individual
        - I_CL: The center line for the individual
        - I_LCL: The lower control limit for the individual
        - MR_UCL: The upper control limit for the moving range
        - MR_CL: The center line for the moving range
        - MR_LCL: The lower control limit for the moving range
        - I_TEST1: The points that violate the alarm rule for the individual
        - MR_TEST1: The points that violate the alarm rule for the moving range
    """
    # Check if df is a pandas DataFrame
    if not isinstance(original_df, pd.DataFrame):
        raise TypeError('The data must be a pandas DataFrame.')

    # Check if the col_name exists in the DataFrame
    if col_name not in original_df.columns:
        raise ValueError('The column name does not exist in the DataFrame.')

    d2 = 1.128
    D4 = 3.267

    # Create a copy of the original DataFrame
    df = original_df.copy()
    
    # Calculate the moving range
    df['MR'] = df[col_name].diff().abs()
    # Create columns for the upper and lower control limits
    df['I_UCL'] = df[col_name].mean() + (3*df['MR'].mean()/d2)
    df['I_CL'] = df[col_name].mean()
    df['I_LCL'] = df[col_name].mean() - (3*df['MR'].mean()/d2)

    # Define columns for the Western Electric alarm rules
    df['I_TEST1'] = np.where((df[col_name] > df['I_UCL']) | (df[col_name] < df['I_LCL']), df[col_name], np.nan)

    # Plot the I and MR charts
    fig, ax = plt.subplots(2, 1)
    fig.suptitle(('I-MR charts of %s' % col_name))
    ax[0].plot(df[col_name], color='mediumblue', linestyle='--', marker='o')
    ax[0].plot(df['I_UCL'], color='firebrick', linewidth=1)
    ax[0].plot(df['I_CL'], color='g', linewidth=1)
    ax[0].plot(df['I_LCL'], color='firebrick', linewidth=1)
    ax[0].set_ylabel('Individual Value')
    ax[1].set_xlabel('Sample number')
    # add the values of the control limits on the right side of the plot
    ax[0].text(len(df)+.5, df['I_UCL'].iloc[0], 'UCL = {:.2f}'.format(df['I_UCL'].iloc[0]), verticalalignment='center')
    ax[0].text(len(df)+.5, df['I_CL'].iloc[0], 'CL = {:.2f}'.format(df['I_CL'].iloc[0]), verticalalignment='center')
    ax[0].text(len(df)+.5, df['I_LCL'].iloc[0], 'LCL = {:.2f}'.format(df['I_LCL'].iloc[0]), verticalalignment='center')
    # highlight the points that violate the alarm rules
    ax[0].plot(df['I_TEST1'], linestyle='none', marker='s', color='firebrick', markersize=10)

    # Create columns for the upper and lower control limits
    df['MR_UCL'] = D4 * df['MR'].mean()
    df['MR_CL'] = df['MR'].mean()
    df['MR_LCL'] = 0
    # Define columns for the Western Electric alarm rules
    df['MR_TEST1'] = np.where((df['MR'] > df['MR_UCL']) | (df['MR'] < df['MR_LCL']), df['MR'], np.nan)

    ax[1].plot(df['MR'], color='mediumblue', linestyle='--', marker='o')
    ax[1].plot(df['MR_UCL'], color='firebrick', linewidth=1)
    ax[1].plot(df['MR_CL'], color='g', linewidth=1)
    ax[1].plot(df['MR_LCL'], color='firebrick', linewidth=1)
    ax[1].set_ylabel('Moving Range')
    ax[1].set_xlabel('Sample number')
    # add the values of the control limits on the right side of the plot
    ax[1].text(len(df)+.5, df['MR_UCL'].iloc[0], 'UCL = {:.2f}'.format(df['MR_UCL'].iloc[0]), verticalalignment='center')
    ax[1].text(len(df)+.5, df['MR_CL'].iloc[0], 'CL = {:.2f}'.format(df['MR_CL'].iloc[0]), verticalalignment='center')
    ax[1].text(len(df)+.5, df['MR_LCL'].iloc[0], 'LCL = {:.2f}'.format(df['MR_LCL'].iloc[0]), verticalalignment='center')
    # highlight the points that violate the alarm rules
    ax[1].plot(df['MR_TEST1'], linestyle='none', marker='s', color='firebrick', markersize=10)


    plt.show()

    return df

import scipy.integrate as spi
import numpy as np

import numpy as np
import scipy.integrate as spi

class constants:
    @staticmethod
    def getd2(n=None):
        if n is None or n < 2 or abs(n - round(n)) != 0:
            raise ValueError("Invalid sample size ({})".format(n))
        def f(x):
            return stats.studentized_range.sf(x, n, np.inf)
        d2, _ = spi.quad(f, 0, np.inf)
        if _ > 0.001:
            print("The absolute error after numerical integration is greater than 0.001")
        return d2

def summary(results):
    """Prints a summary of the regression results.

    Parameters
    ----------
    results : RegressionResults object
        The results of a regression model.

    Returns
    -------
    None
    """

    # Set the precision of the output
    np.set_printoptions(precision=4, suppress=True)
    pd.options.display.precision = 4

    # Extract information from the result object
    terms = results.model.exog_names
    coefficients = results.params
    std_errors = results.bse
    t_values = results.tvalues
    p_values = results.pvalues
    #r_squared = results.rsquared
    #adjusted_r_squared = results.rsquared_adj

    # Print the regression equation
    print("REGRESSION EQUATION")
    print("-------------------")
    equation = ("%s = " % results.model.endog_names)
    for i in range(len(coefficients)):
        if results.model.exog_names[i] == 'Intercept':
            equation += "%.3f" % coefficients[i]
        else:
            if coefficients[i] > 0:
                equation += " + %.3f %s" % (coefficients[i], terms[i])
            else:
                equation += " %.3f %s" % (coefficients[i], terms[i])
    print(equation)

    # Print the information in a similar format to Minitab
    print("\nCOEFFICIENTS")
    print("------------")
    # make a dataframe to store the results
    
    df_coefficients = pd.DataFrame({'Term': terms, 'Coef': coefficients, 'SE Coef': std_errors, 'T-Value': t_values, 'P-Value': p_values})
    df_coefficients.style.format({'Coef': '{:.3f}', 'SE Coef': '{:.3f}', 'T-Value': '{:.3f}', 'P-Value': '{:.3f}'})
    print(df_coefficients.to_string(index=False))

    # Print the R-squared and adjusted R-squared
    print("\nMODEL SUMMARY")
    print("-------------")
    # compute the standard deviation of the distance between the data values and the fitted values
    S = np.std(results.resid, ddof=len(terms))
    # make a dataframe to store the results
    df_model_summary = pd.DataFrame({'S': [S], 'R-sq': [results.rsquared], 'R-sq(adj)': [results.rsquared_adj]})
    print(df_model_summary.to_string(index=False))

    # Print the ANOVA table
    print("\nANALYSIS OF VARIANCE")
    print("---------------------")
    # make a dataframe with the column names and no data
    df_anova = pd.DataFrame(columns=['Source', 'DF', 'Adj SS', 'Adj MS', 'F-Value', 'P-Value'])
    # add the rows of data
    df_anova.loc[0] = ['Regression', results.df_model, results.mse_model * results.df_model, results.mse_model, results.fvalue, results.f_pvalue]
    jj = 1
    for term in terms:
        if term != 'Intercept':
            # perform the f-test for the term
            f_test = results.f_test(term + '= 0')
            df_anova.loc[jj] = [term, f_test.df_num, f_test.fvalue * results.mse_resid * f_test.df_num, f_test.fvalue * results.mse_resid, f_test.fvalue, f_test.pvalue]
            jj += 1

    df_anova.loc[jj] = ['Error', results.df_resid, results.mse_resid * results.df_resid, results.mse_resid, np.nan, np.nan]

    '''
    # Lack-of-fit
    # compute the number of levels in the independent variables 
    n_levels = results.df_resid
    for term in terms:
        if term != 'Intercept':
            n_levels = np.minimum(n_levels, len(data[term].unique())

    if n_levels < results.df_resid:
        dof_lof = n_levels - len(terms)
        dof_pe = results.df_resid - n_levels
        # compute the SSE for the pure error term
        for 


        df_anova.loc[jj + 1] = ['Lack-of-fit', n_levels - len(terms), np.nan, np.nan, np.nan, np.nan]
    '''

    df_anova.loc[jj + 1] = ['Total', results.df_model + results.df_resid, results.mse_total * (results.df_model + results.df_resid), np.nan, np.nan, np.nan]

    print(df_anova.to_string(index=False))

    return


def ARIMAsummary(results):

    # Set the precision of the output
    np.set_printoptions(precision=4, suppress=True)
    pd.options.display.precision = 4

    # Extract information from the result object
    terms = results.param_names
    coefficients = results.params
    std_errors = results.bse
    t_values = results.tvalues
    p_values = results.pvalues
    n_coefficients = len(coefficients) - 1 #because models givers an additional information on sigma^2 in the list of coefficients
    
    
    # get the order of the model
    n_model = results.nobs
    ar_order = results.model.order[0]
    ma_order = results.model.order[2]
    diff_order = results.model.order[1]
    order_model = results.model.order
    order_model_flag = sum(order_model) > 0
    
    #get seasonal order vector
    so_model = results.model.seasonal_order
    DIFF_seasonal_order = so_model[1]
    seasonal_model_flag = so_model[3] > 0

    #Model's degrees of freedom
    df_model = (results.nobs - diff_order - DIFF_seasonal_order) - len(results.params) - 1 #degrees of freedom for the model: (n - d - D) - estimated parameters(p, q, P, Q, constant term)

    print("---------------------")
    print("ARIMA MODEL RESULTS")
    print("---------------------")

    if order_model_flag:
        print(f"ARIMA model order: p={ar_order}, d={diff_order}, q={ma_order}")
    if seasonal_model_flag:
        print(f"Seasonal ARIMA model fit with period {so_model[3]} and order: P={so_model[0]}, D={so_model[1]}, Q={so_model[2]}")

    
    # Print the information in a similar format to Minitab
    print("\nFINAL ESTIMATES OF PARAMETERS")
    print("-------------------------------")
    # make a dataframe to store the results
    
    df_coefficients = pd.DataFrame({'Term': terms[0:n_coefficients], 'Coef': coefficients[0:n_coefficients], 'SE Coef': std_errors[0:n_coefficients], 'T-Value': t_values[0:n_coefficients], 'P-Value': p_values[0:n_coefficients]})
    df_coefficients.style.format({'Coef': '{:.3f}', 'SE Coef': '{:.3f}', 'T-Value': '{:.3f}', 'P-Value': '{:.3f}'})
    print(df_coefficients.to_string(index=False))


    # Print the ANOVA table
    print("\nRESIDUAL SUM OF SQUARES")
    print("-------------------------")
    # make a dataframe with the column names and no data
    df_rss = pd.DataFrame(columns=['DF', 'SS', 'MS'])
    # add the rows of data
    df_model = len(results.fittedvalues) - results.df_model + 1 #degrees of freedom= n-(d+p+q+1 [if constant term is included]])
    df_rss.loc[0] = [df_model, results.sse, results.mse]
    print(df_rss.to_string(index=False))


    # Print the information in a similar format to Minitab for LBQ test
    print("\nLjung-Box Chi-Square Statistics")
    print("----------------------------------")
    if len(results.resid) > 48:
        lagvalues = np.array([12, 24, 36, 48])
    elif len(results.resid) > 36:
        lagvalues = np.array([12, 24, 36])
    elif len(results.resid) > 24:
        lagvalues = np.array([12, 24])
    elif len(results.resid) > 12:
        lagvalues = np.array([12])
    else:
        lagvalues = int(np.sqrt(len(results.resid)))
    LBQ=acorr_ljungbox(results.resid, lags=lagvalues, return_df=True, boxpierce=True)

    df_LBtest = pd.DataFrame({'Lag': lagvalues, 'Chi-Square': LBQ.lb_stat, 'P-Value': LBQ.lb_pvalue})
    df_LBtest.style.format({'Lag': '{:.3f}', 'Chi-Square': '{:.3f}', 'P-Value': '{:.3f}'})
    print(df_LBtest.to_string(index=False))
   
    return


# create a class called StepwiseRegression that performs stepwise regression when fitting a model
class StepwiseRegression:
    
    """Performs stepwise regression.

    Parameters
    ----------
    
    y : array-like
        The dependent variable.
    X : array-like
        The independent variables.
    
    add_constant : bool, optional
        Whether to add a constant to the model. The default is True.
    direction : string, optional
        The direction of stepwise regression. The default is 'both'.
    alpha_to_enter : float, optional
        The alpha level to enter a variable in the forward step. The default is 0.15.
    alpha_to_remove : float, optional
        The alpha level to remove a variable in the backward step. The default is 0.15.
    max_iterations : int, optional
        The maximum number of iterations. The default is 100.

    Returns
    -------
    model_fit : RegressionResults object
        The results of a regression model.

    """

    # initialize the class
    def __init__(self, add_constant = True, direction = 'both', alpha_to_enter = 0.15, alpha_to_remove = 0.15, max_iterations = 100):
        self.add_constant = add_constant
        self.direction = direction
        self.alpha_to_enter = alpha_to_enter
        self.alpha_to_remove = alpha_to_remove
        self.max_iterations = max_iterations
        self.break_loop = False
        self.model_fit = None

    # define a function to fit the model
    def fit(self, y, X):
        self.X = X
        self.y = y
        self.variables_to_include = []
        
        # fit the initial models with one independent variable at a time
        print('Stepwise Regression')
        print('\n######################################')
        k = 1
        print('### Step %d' % k)
        print('-------------------')
        self.forward_selection()

        # check if self.variables_to_include is empty
        if len(self.variables_to_include) == 0:
            raise ValueError('All predictors have p-values greater than the alpha_to_enter level. No model was selected.')
        
        while self.break_loop == False:
            k += 1
            print('\n######################################')
            print('### Step %d' % k)
            print('-------------------')
            if self.direction == 'both':
                self.forward_selection()
                print('-------------------')
                if self.break_loop == False:
                    self.backward_elimination()
            else:
                raise ValueError('The direction must be either "both", "forward", or "backward".')
            
            if k == self.max_iterations:
                self.break_loop = True
                print('Maximum number of iterations reached.')

        return self        
    
    def forward_selection(self):

        print('Forward Selection')

        selected_pvalue = self.alpha_to_enter
        if len(self.variables_to_include) == 0:
            original_variables = []
        else:
            original_variables = self.variables_to_include
        
        number_of_variables = len(self.variables_to_include)

        if number_of_variables == self.X.shape[1]:
            self.break_loop = True
            print('All predictors have been included in the model. Exiting stepwise.')
            return self

        # fit the model with the selected variables and add one of the remaining variables at a time
        for i in range(self.X.shape[1]):

            if i not in self.variables_to_include:
                # create a new list called testing_variables that includes the original variables and the new variable
                testing_variables = original_variables.copy()
                testing_variables.append(i)

                X_test = self.X.iloc[:, testing_variables]
                
                if self.add_constant:
                    X_test = sm.add_constant(X_test)

                model_fit = sm.OLS(self.y, X_test).fit()
                
                # if the p-value of the new variable is less than the alpha_to_enter level, 
                # add the variable to the list of variables to include
                if model_fit.pvalues[-1] < self.alpha_to_enter and model_fit.pvalues[-1] < selected_pvalue:
                    selected_pvalue = model_fit.pvalues[-1]
                    self.variables_to_include = testing_variables
                    self.model_fit = model_fit

        if len(self.variables_to_include) == number_of_variables:
            self.break_loop = True
            print('\nNo predictor added. Exiting stepwise.')
        else:
            # print(self.model_fit.summary())
            self.SWsummary()

        return self
        
    
    def backward_elimination(self):
        
        print('Backward Elimination')
        
        original_variables = self.variables_to_include

        # sort the pvalues in descending order and remove the variable with pvalue > alpha_to_remove
        if self.add_constant:
            sorted_pvalues = self.model_fit.pvalues[1:].sort_values(ascending = False)
        else:
            sorted_pvalues = self.model_fit.pvalues.sort_values(ascending = False)

        testing_variables = original_variables.copy()

        for i in range(len(sorted_pvalues)):
            if sorted_pvalues[i] > self.alpha_to_remove:
                variable_to_remove = sorted_pvalues.index[i]
                testing_variables.remove(self.X.columns.get_loc(variable_to_remove))
            else:
                break
        
        if len(testing_variables) == len(original_variables):
            print('\nNo predictor removed.')
            return(self)
        
        X_test = self.X.iloc[:, testing_variables]

        if self.add_constant:
            X_test = sm.add_constant(X_test)

        self.model_fit = sm.OLS(self.y, X_test).fit()
        self.SWsummary()

        return self
    
    def SWsummary(self):
        # Extract information from the result object
        results = self.model_fit
        terms = results.model.exog_names
        coefficients = results.params
        p_values = results.pvalues
        #r_squared = results.rsquared
        #adjusted_r_squared = results.rsquared_adj

        # Print the information in a similar format to Minitab
        print("\nCOEFFICIENTS")
        print("------------")
        # make a dataframe to store the results
        
        df_coefficients = pd.DataFrame({'Term': terms, 'Coef': coefficients, 'P-Value': p_values})
        print(df_coefficients.to_string(index=False))

        # Print the R-squared and adjusted R-squared
        print("\nMODEL SUMMARY")
        print("-------------")
        # compute the standard deviation of the distance between the data values and the fitted values
        S = np.std(results.resid, ddof=len(terms))
        # make a dataframe to store the results
        df_model_summary = pd.DataFrame({'S': [S], 'R-sq': [results.rsquared], 'R-sq(adj)': [results.rsquared_adj]})
        print(df_model_summary.to_string(index=False))

        return        


# #%% create a stepwise regression object
# stepwise = StepwiseRegression(direction = 'both', alpha_to_enter = 0.15, alpha_to_remove = 0.15)

# # create a dataframe with the independent variables
# np.random.seed(1)
# x1 = np.arange(100)
# x2 = np.exp(np.arange(100)*0.03)
# x3 = 2 * x1
# x4 = np.random.normal(9, 1, 100)
# X = pd.DataFrame({'x1': x1, 'x2': x2, 'x3': x3, 'x4': x4})
# y = pd.DataFrame({'y': 0.5*x1 + 5*x2 + np.random.normal(0, 10, 100)})

# plt.plot(x1, y, 'o')
# plt.xlabel('x1')
# plt.ylabel('y')
# plt.show()

# plt.plot(x2, y, 'o')
# plt.xlabel('x2')
# plt.ylabel('y')
# plt.show()

# #%%
# model = stepwise.fit(X, y)


# # %%
# #Import the necessary libraries
# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd
# from scipy import stats

# #Import the dataset
# data = pd.read_csv('ESE3_es5_dataset.csv')

# # Now add a 'year' column to the dataset with values from 1900 to 1997
# data['year'] = np.arange(1900, 1998)
# print(data.head())

# # Add a column to the dataset with the lagged values
# data['Ex5_lag1'] = data['Ex5'].shift(1)
# # Add a lag4 term to the dataframe
# data['Ex5_lag4'] = data['Ex5'].shift(4)

# # and split the dataset into regressors and target
# X = data.iloc[4:, 1:]
# y = data.iloc[4:, 0]

# stepwise = StepwiseRegression(add_constant = True, direction = 'both', alpha_to_enter = 0.15, alpha_to_remove = 0.15)

# # Fit the model
# model = stepwise.fit(y, X)

# print(summary(model.model_fit))

# # reorder the model coefficients to match the order of the variables in the original dataset
# X = sm.add_constant(X)
# model.model_fit.params = model.model_fit.params.reindex(X.columns)

# print(summary(model.model_fit))
