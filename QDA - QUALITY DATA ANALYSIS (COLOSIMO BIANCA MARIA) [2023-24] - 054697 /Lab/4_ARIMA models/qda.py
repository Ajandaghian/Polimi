#%% version 1.5
__version__ = 'qda version 1.5'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.diagnostic import acorr_ljungbox
import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA as arimafromlib


class ControlCharts:
    @staticmethod
    def IMR(original_df, col_name, K = 3, subset_size = None, run_rules = False, plotit = True):
        """Implements the Individual Moving Range (IMR) chart.
        Parameters
        ----------
        original_df : pandas.DataFrame
            A DataFrame containing the data to be plotted.
        col_name : str
            The name of the column to be used for the IMR control chart.
        K : int, optional
            The number of standard deviations. Default is 3.
        subset_size : int, optional
            The number of rows to be used for the IMR chart. Default is None and all rows are used.

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

        # get the size of the DataFrame
        n = 2
        d2 = constants.getd2(n)
        D4 = constants.getD4(n, K)

        if subset_size is None:
            subset_size = len(original_df)
        elif subset_size > len(original_df):
            raise ValueError('The subset size must be less than the number of rows in the DataFrame.')

        # Create a copy of the original DataFrame
        df = original_df.copy()
        
        # Calculate the moving range
        df['MR'] = df[col_name].diff().abs()
        # Create columns for the upper and lower control limits
        df['I_UCL'] = df[col_name].iloc[:subset_size].mean() + (K*df['MR'].iloc[:subset_size].mean()/d2)
        df['I_CL'] = df[col_name].iloc[:subset_size].mean()
        df['I_LCL'] = df[col_name].iloc[:subset_size].mean() - (K*df['MR'].iloc[:subset_size].mean()/d2)
        # Define columns for the Western Electric alarm rules
        df['I_TEST1'] = np.where((df[col_name] > df['I_UCL']) | (df[col_name] < df['I_LCL']), df[col_name], np.nan)

        # Create columns for the upper and lower control limits
        df['MR_UCL'] = D4 * df['MR'].iloc[:subset_size].mean()
        df['MR_CL'] = df['MR'].iloc[:subset_size].mean()
        df['MR_LCL'] = 0
        # Define columns for the Western Electric alarm rules
        df['MR_TEST1'] = np.where((df['MR'] > df['MR_UCL']) | (df['MR'] < df['MR_LCL']), df['MR'], np.nan)

        if plotit == True:
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

            if subset_size < len(original_df):
                ax[0].axvline(x=subset_size-.5, color='k', linestyle='--')
                ax[1].axvline(x=subset_size-.5, color='k', linestyle='--')

        

        if run_rules == True:

            # Create a new dataframe for the Run Rules
            df_rr = pd.DataFrame()

            # add columns for the RR calculation
            # RR2 count the number of consecutive points on the same side of the center line
            K_2 = 9 
            df['I_TEST2'] = np.nan
            df_rr['I_TEST2_1'] = df[col_name] - df['I_CL']
            df_rr['I_TEST2_2'] = 0
            df_rr['I_TEST2_2'][0] = 1

            # RR3 count the number of consecutive increasing or decreasing points
            K_3 = 6
            df_IMR['I_TEST3'] = np.nan
            df_rr['I_TEST3_1'] = df_IMR[col_name].diff()
            df_rr['I_TEST3_2'] = 0
            df_rr['I_TEST3_2'][0] = 1

            # RR4 count the number of consecutive points alternating up and down around the center line
            K_4 = 14
            df['I_TEST4'] = np.nan
            df_rr['I_TEST4_1'] = 0
            df_rr['I_TEST4_1'][0] = 1

            # RR5 K out of K+1 points > 2 standard deviations from the center line (same side)
            sigma_I = df['MR'].mean()/d2
            
            K_5 = 2
            df['I_TEST5'] = np.nan
            df_rr['I_TEST5_1'] = (df[col_name] - (df['I_CL'] + 2 * sigma_I)) > 0
            df_rr['I_TEST5_2'] = (df[col_name] - (df['I_CL'] - 2 * sigma_I)) < 0
            # create a moving average of the previous K_5 points
            df_rr['I_TEST5_3'] = df_rr['I_TEST5_1'].rolling(K_5+1).mean()
            df_rr['I_TEST5_4'] = df_rr['I_TEST5_2'].rolling(K_5+1).mean()

            # RR6 K out of K+1 points > 1 standard deviations from the center line (same side)
            K_6 = 4
            df['I_TEST6'] = np.nan
            df_rr['I_TEST6_1'] = (df[col_name] - (df['I_CL'] + 1 * sigma_I)) > 0
            df_rr['I_TEST6_2'] = (df[col_name] - (df['I_CL'] - 1 * sigma_I)) < 0
            # create a moving average of the previous K_5 points
            df_rr['I_TEST6_3'] = df_rr['I_TEST5_1'].rolling(K_6+1).mean()
            df_rr['I_TEST6_4'] = df_rr['I_TEST5_2'].rolling(K_6+1).mean()

            # RR7 K out of K+1 points within 1 standard deviations from the center line (either side)
            K_7 = 15
            df['I_TEST7'] = np.nan
            # create a column that returns 1 if the point is within 1 standard deviation of the center line 
            df_rr['I_TEST7_1'] = np.abs(df[col_name] - df['I_CL']) < sigma_I
            # create a moving average of the previous K_7 points
            df_rr['I_TEST7_2'] = df_rr['I_TEST7_1'].rolling(K_7+1).mean()

            # RR8 K out of K+1 points > 1 standard deviations from the center line (either side)
            K_8 = 8
            df['I_TEST8'] = np.nan
            # create a column that returns 1 if the point is > 1 standard deviation of the center line
            df_rr['I_TEST8_1'] = np.abs(df[col_name] - df['I_CL']) > sigma_I
            # create a moving average of the previous K_8 points
            df_rr['I_TEST8_2'] = df_rr['I_TEST8_1'].rolling(K_8+1).mean()



            for i in range(1, len(df_IMR)):
                # Rule 2
                if df_rr['I_TEST2_1'][i] * df_rr['I_TEST2_1'][i-1] > 0:
                    df_rr['I_TEST2_2'][i] = df_rr['I_TEST2_2'][i-1] + 1
                else:
                    df_rr['I_TEST2_2'][i] = 1

                if df_rr['I_TEST2_2'][i] >= K_2:
                    df['I_TEST2'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '2', color='k', fontsize=12)

                # Rule 3
                if df_rr['I_TEST3_1'][i] * df_rr['I_TEST3_1'][i-1] > 0:
                    df_rr['I_TEST3_2'][i] = df_rr['I_TEST3_2'][i-1] + 1
                else:
                    df_rr['I_TEST3_2'][i] = 1

                if df_rr['I_TEST3_2'][i] >= K_3:
                    df['I_TEST3'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '3', color='k', fontsize=12)

                # Rule 4
                if df_rr['I_TEST3_1'][i] * df_rr['I_TEST3_1'][i-1] < 0:
                    df_rr['I_TEST4_1'][i] = df_rr['I_TEST4_1'][i-1] + 1
                else:
                    df_rr['I_TEST4_1'][i] = 1

                if df_rr['I_TEST4_1'][i] >= K_4:
                    df['I_TEST4'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '4', color='k', fontsize=12)

                # Rule 5
                if df_rr['I_TEST5_3'][i] >= K_5/(K_5+1):
                    df['I_TEST5'][i] = df[col_name][i]
                elif df_rr['I_TEST5_4'][i] >= K_5/(K_5+1):
                    df['I_TEST5'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '5', color='k', fontsize=12)
                
                # Rule 6
                if df_rr['I_TEST6_3'][i] >= K_6/(K_6+1):
                    df['I_TEST6'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '6', color='k', fontsize=12)
                elif df_rr['I_TEST6_4'][i] >= K_6/(K_6+1):
                    df['I_TEST6'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '6', color='k', fontsize=12)

                # Rule 7
                if df_rr['I_TEST7_2'][i] >= K_7/(K_7+1):
                    df['I_TEST7'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '7', color='k', fontsize=12)

                # Rule 8
                if df_rr['I_TEST8_2'][i] >= K_8/(K_8+1):
                    df['I_TEST8'][i] = df[col_name][i]
                    ax[0].text(i, df[col_name][i], '7', color='k', fontsize=12)
            
            # Plot the alarm rules
            ax[0].plot(df['I_TEST2'], linestyle='none', marker='X', color='orange', markersize=10)
            ax[0].plot(df['I_TEST3'], linestyle='none', marker='X', color='orange', markersize=10)
            ax[0].plot(df['I_TEST4'], linestyle='none', marker='X', color='orange', markersize=10)
            ax[0].plot(df['I_TEST5'], linestyle='none', marker='X', color='orange', markersize=10)
            ax[0].plot(df['I_TEST6'], linestyle='none', marker='X', color='orange', markersize=10)
            ax[0].plot(df['I_TEST7'], linestyle='none', marker='X', color='orange', markersize=10)
            ax[0].plot(df['I_TEST8'], linestyle='none', marker='X', color='orange', markersize=10)
            
        plt.show()

        return df
    
    @staticmethod
    def XbarR(original_df, K = 3, subset_size = None, plotit = True):
        '''
        This function plots the Xbar-R charts of a DataFrame 
        and returns the DataFrame with the control limits and alarm rules.

        Parameters
        ----------
        original_df : DataFrame
            The DataFrame that contains the data.
        K : int, optional
            The number of standard deviations. The default is 3.
        subset_size : int, optional
            The number of rows to be used for the IMR chart. Default is None and all rows are used.

        Returns
        -------
        data_XR : DataFrame
            The DataFrame with the control limits and alarm rules.
        '''
        # get the shape of the DataFrame
        m, n = original_df.shape

        if n < 2:
            raise ValueError('The DataFrame must contain at least 2 columns.')

        # Calculate the constants
        A2 = constants.getA2(n, K)
        D4 = constants.getD4(n, K)
        D3 = constants.getD3(n, K)

        if subset_size is None:
            subset_size = len(original_df)
        elif subset_size > len(original_df):
            raise ValueError('The subset size must be less than the number of rows in the DataFrame.')

        # Create a copy of the original DataFrame
        data_XR = original_df.copy()

        # Add a column with the mean of the rows
        data_XR['sample_mean'] = original_df.mean(axis=1)
        # Add a column with the range of the rows
        data_XR['sample_range'] = original_df.max(axis=1) - original_df.min(axis=1)

        Xbar_mean = data_XR['sample_mean'].iloc[:subset_size].mean()
        R_mean = data_XR['sample_range'].iloc[:subset_size].mean()

        # Now we can compute the CL, UCL and LCL for Xbar and R
        data_XR['Xbar_CL'] = Xbar_mean
        data_XR['Xbar_UCL'] = Xbar_mean + A2 * R_mean
        data_XR['Xbar_LCL'] = Xbar_mean - A2 * R_mean

        data_XR['R_CL'] = R_mean
        data_XR['R_UCL'] = D4 * R_mean
        data_XR['R_LCL'] = D3 * R_mean

        # Define columns for the alarms
        data_XR['Xbar_TEST1'] = np.where((data_XR['sample_mean'] > data_XR['Xbar_UCL']) | 
                (data_XR['sample_mean'] < data_XR['Xbar_LCL']), data_XR['sample_mean'], np.nan)
        data_XR['R_TEST1'] = np.where((data_XR['sample_range'] > data_XR['R_UCL']) | 
                (data_XR['sample_range'] < data_XR['R_LCL']), data_XR['sample_range'], np.nan)

        fig, ax = plt.subplots(2, 1)
        fig.suptitle(('Xbar-R charts'))
        ax[0].plot(data_XR['sample_mean'], color='mediumblue', linestyle='--', marker='o')
        ax[0].plot(data_XR['Xbar_UCL'], color='firebrick', linewidth=1)
        ax[0].plot(data_XR['Xbar_CL'], color='g', linewidth=1)
        ax[0].plot(data_XR['Xbar_LCL'], color='firebrick', linewidth=1)
        ax[0].set_ylabel('Sample Mean')
        # add the values of the control limits on the right side of the plot
        ax[0].text(len(data_XR)+.5, data_XR['Xbar_UCL'].iloc[0], 'UCL = {:.3f}'.format(data_XR['Xbar_UCL'].iloc[0]), verticalalignment='center')
        ax[0].text(len(data_XR)+.5, data_XR['Xbar_CL'].iloc[0], 'CL = {:.3f}'.format(data_XR['Xbar_CL'].iloc[0]), verticalalignment='center')
        ax[0].text(len(data_XR)+.5, data_XR['Xbar_LCL'].iloc[0], 'LCL = {:.3f}'.format(data_XR['Xbar_LCL'].iloc[0]), verticalalignment='center')
        # highlight the points that violate the alarm rules
        ax[0].plot(data_XR['Xbar_TEST1'], linestyle='none', marker='s', color='firebrick', markersize=10)

        ax[1].plot(data_XR['sample_range'], color='mediumblue', linestyle='--', marker='o')
        ax[1].plot(data_XR['R_UCL'], color='firebrick', linewidth=1)
        ax[1].plot(data_XR['R_CL'], color='g', linewidth=1)
        ax[1].plot(data_XR['R_LCL'], color='firebrick', linewidth=1)
        ax[1].set_ylabel('Sample Range')
        ax[1].set_xlabel('Sample Number')
        # add the values of the control limits on the right side of the plot
        ax[1].text(len(data_XR)+.5, data_XR['R_UCL'].iloc[0], 'UCL = {:.3f}'.format(data_XR['R_UCL'].iloc[0]), verticalalignment='center')
        ax[1].text(len(data_XR)+.5, data_XR['R_CL'].iloc[0], 'CL = {:.3f}'.format(data_XR['R_CL'].iloc[0]), verticalalignment='center')
        ax[1].text(len(data_XR)+.5, data_XR['R_LCL'].iloc[0], 'LCL = {:.3f}'.format(data_XR['R_LCL'].iloc[0]), verticalalignment='center')
        # highlight the points that violate the alarm rules
        ax[1].plot(data_XR['R_TEST1'], linestyle='none', marker='s', color='firebrick', markersize=10)

        if subset_size < len(original_df):
            ax[0].axvline(x=subset_size-.5, color='k', linestyle='--')
            ax[1].axvline(x=subset_size-.5, color='k', linestyle='--')

        plt.show()

        return data_XR
    
    @staticmethod
    def XbarS(original_df, K = 3, sigma = None, subset_size = None, plotit = True):
        '''
        This function plots the Xbar-S charts of a DataFrame 
        and returns the DataFrame with the control limits and alarm rules.

        Parameters
        ----------
        original_df : DataFrame
            The DataFrame that contains the data.
        K : int, optional
            The number of standard deviations. The default is 3.
        sigma : float, optional
            Input the standard deviation of the population. Otherwise, the standard deviation of the sample will be used.
        subset_size : int, optional
            The number of rows to be used for the IMR chart. Default is None and all rows are used.

        Returns
        -------
        data_XS : DataFrame
            The DataFrame with the control limits and alarm rules.
        '''
        # get the shape of the DataFrame
        m, n = original_df.shape

        if n < 2:
            raise ValueError('The DataFrame must contain at least 2 columns.')
        
        # Choose the data subset to be used for the control limits
        if subset_size is None:
            subset_size = len(original_df)
        elif subset_size > len(original_df):
            raise ValueError('The subset size must be less than the number of rows in the DataFrame.')

        # Calculate the constants
        if sigma is None:
            A3 = K * 1 / (constants.getc4(n) * np.sqrt(n))
            B3 = np.maximum(1 - K * (np.sqrt(1-constants.getc4(n)**2)) / (constants.getc4(n)), 0)
            B4 = 1 + K * (np.sqrt(1-constants.getc4(n)**2)) / (constants.getc4(n))
        else:
            A3 = K * 1 / (constants.getc4(n) * np.sqrt(n))
            B5 = np.maximum(constants.getc4(n) - K * np.sqrt(1-constants.getc4(n)**2), 0)
            B6 = constants.getc4(n) + K * np.sqrt(1-constants.getc4(n)**2)

        # Create a copy of the original DataFrame
        data_XS = original_df.copy()

        # Add a column with the mean of the rows
        data_XS['sample_mean'] = original_df.mean(axis=1)
        # Add a column with the range of the rows
        data_XS['sample_std'] = original_df.std(axis=1)

        Xbar_mean = data_XS['sample_mean'].iloc[:subset_size].mean()
        S_mean = data_XS['sample_std'].iloc[:subset_size].mean()

        # Now we can compute the CL, UCL and LCL for Xbar and S
        data_XS['Xbar_CL'] = Xbar_mean
        data_XS['Xbar_UCL'] = Xbar_mean + A3 * S_mean
        data_XS['Xbar_LCL'] = Xbar_mean - A3 * S_mean

        if sigma is None:
            data_XS['S_CL'] = S_mean
            data_XS['S_UCL'] = B4 * S_mean
            data_XS['S_LCL'] = B3 * S_mean
        else:
            data_XS['S_CL'] = sigma
            data_XS['S_UCL'] = B6 * sigma
            data_XS['S_LCL'] = B5 * sigma

        # Define columns for the alarms
        data_XS['Xbar_TEST1'] = np.where((data_XS['sample_mean'] > data_XS['Xbar_UCL']) | 
                        (data_XS['sample_mean'] < data_XS['Xbar_LCL']), data_XS['sample_mean'], np.nan)
        data_XS['S_TEST1'] = np.where((data_XS['sample_std'] > data_XS['S_UCL']) | 
                        (data_XS['sample_std'] < data_XS['S_LCL']), data_XS['sample_std'], np.nan)

        fig, ax = plt.subplots(2, 1)
        fig.suptitle(('Xbar-S charts'))
        ax[0].plot(data_XS['sample_mean'], color='mediumblue', linestyle='--', marker='o')
        ax[0].plot(data_XS['Xbar_UCL'], color='firebrick', linewidth=1)
        ax[0].plot(data_XS['Xbar_CL'], color='g', linewidth=1)
        ax[0].plot(data_XS['Xbar_LCL'], color='firebrick', linewidth=1)
        ax[0].set_ylabel('Sample Mean')
        # add the values of the control limits on the right side of the plot
        ax[0].text(len(data_XS)+.5, data_XS['Xbar_UCL'].iloc[0], 'UCL = {:.3f}'.format(data_XS['Xbar_UCL'].iloc[0]), verticalalignment='center')
        ax[0].text(len(data_XS)+.5, data_XS['Xbar_CL'].iloc[0], 'CL = {:.3f}'.format(data_XS['Xbar_CL'].iloc[0]), verticalalignment='center')
        ax[0].text(len(data_XS)+.5, data_XS['Xbar_LCL'].iloc[0], 'LCL = {:.3f}'.format(data_XS['Xbar_LCL'].iloc[0]), verticalalignment='center')
        # highlight the points that violate the alarm rules
        ax[0].plot(data_XS['Xbar_TEST1'], linestyle='none', marker='s', color='firebrick', markersize=10)

        ax[1].plot(data_XS['sample_std'], color='mediumblue', linestyle='--', marker='o')
        ax[1].plot(data_XS['S_UCL'], color='firebrick', linewidth=1)
        ax[1].plot(data_XS['S_CL'], color='g', linewidth=1)
        ax[1].plot(data_XS['S_LCL'], color='firebrick', linewidth=1)
        ax[1].set_ylabel('Sample StDev')
        ax[1].set_xlabel('Sample Number')
        # add the values of the control limits on the right side of the plot
        ax[1].text(len(data_XS)+.5, data_XS['S_UCL'].iloc[0], 'UCL = {:.3f}'.format(data_XS['S_UCL'].iloc[0]), verticalalignment='center')
        ax[1].text(len(data_XS)+.5, data_XS['S_CL'].iloc[0], 'CL = {:.3f}'.format(data_XS['S_CL'].iloc[0]), verticalalignment='center')
        ax[1].text(len(data_XS)+.5, data_XS['S_LCL'].iloc[0], 'LCL = {:.3f}'.format(data_XS['S_LCL'].iloc[0]), verticalalignment='center')
        # highlight the points that violate the alarm rules
        ax[1].plot(data_XS['S_TEST1'], linestyle='none', marker='s', color='firebrick', markersize=10)

        if subset_size < len(original_df):
            ax[0].axvline(x=subset_size-.5, color='k', linestyle='--')
            ax[1].axvline(x=subset_size-.5, color='k', linestyle='--')

        
        plt.show()

        return data_XS


import numpy as np
import scipy.integrate as spi
import scipy.special as sps

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

    @staticmethod
    def getd3(n=None):
        if n is None or n < 2 or abs(n - round(n)) != 0:
            raise ValueError("Invalid sample size")
        def f(x):
            return x * stats.studentized_range.sf(x, n, np.inf)
        d3, _ = spi.quad(f, 0, np.inf)
        if _ > 0.001:
            print("The absolute error after numerical integration is greater than 0.001")
        d3 = 2 * d3
        this_d2 = constants.getd2(n)
        d3 = np.sqrt(d3 - this_d2**2)
        return d3

    @staticmethod
    def getc4(n=None):
        if n is None or n < 2 or abs(n - round(n)) != 0:
            raise ValueError("Invalid sample size")
        c4 = np.sqrt(2 / (n-1)) * (sps.gamma(n/2) / sps.gamma((n-1)/2))
        return c4

    @staticmethod
    def getA2(n=None, K = 3):
        if n is None or n < 2 or abs(n - round(n)) != 0:
            raise ValueError("Invalid sample size")
        A2 = K / (constants.getd2(n) * np.sqrt(n))
        return A2

    @staticmethod
    def getD3(n=None, K = 3):
        if n is None or n < 2 or abs(n - round(n)) != 0:
            raise ValueError("Invalid sample size")
        D3 = np.maximum(0, 1 - K * constants.getd3(n) / constants.getd2(n))
        return D3

    @staticmethod
    def getD4(n=None, K = 3):
        if n is None or n < 2 or abs(n - round(n)) != 0:
            raise ValueError("Invalid sample size")
        D4 = 1 + K * constants.getd3(n) / constants.getd2(n)
        return D4

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

    """Prints a summary of the ARIMA results.

    Parameters
    ----------
    results : ARIMA object
        The results of an ARIMA.

    Returns
    -------
    None
    """

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
    max_order=np.max(results.model.order)
    
    #get seasonal order vector
    so_model = results.model.seasonal_order
    DIFF_seasonal_order = so_model[1]
    seasonal_model_flag = so_model[3] > 0
    

    #Model's degrees of freedom
    df_model = (results.nobs - diff_order - DIFF_seasonal_order) - (len(results.params) - 1) #degrees of freedom for the model: (n - d - D) - estimated parameters(p, q, P, Q, constant term)

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
    SSE = np.sum(results.resid[max_order:]**2)

    df_rss.loc[0] = [df_model, SSE, SSE/df_model]
    print(df_rss.to_string(index=False))


    # Print the information in a similar format to Minitab for LBQ test
    print("\nLjung-Box Chi-Square Statistics")
    print("----------------------------------")
    if len(results.resid[max_order:]) > 48:
        lagvalues = np.array([12, 24, 36, 48])
    elif len(results.resid[max_order:]) > 36:
        lagvalues = np.array([12, 24, 36])
    elif len(results.resid[max_order:]) > 24:
        lagvalues = np.array([12, 24])
    elif len(results.resid[max_order:]) > 12:
        lagvalues = np.array([12])
    else:
        lagvalues = int(np.sqrt(len(results.resid[max_order:])))
    LBQ=acorr_ljungbox(results.resid[max_order:], lags=lagvalues, boxpierce=True)

    df_LBtest = pd.DataFrame({'Lag': lagvalues, 'Chi-Square': LBQ.lb_stat, 'P-Value': LBQ.lb_pvalue})
    df_LBtest.style.format({'Lag': '{:.3f}', 'Chi-Square test': '{:.3f}', 'P-Value': '{:.3f}'})
    print(df_LBtest.to_string(index=False))
   
    return


def ARIMA(x, order, add_constant):
    """Fits an ARIMA model.

    Parameters
    ----------
    x : data object
    
    order : tuple
        The order of the ARIMA model as (p, d, q)

    add_constant : bool
        True if the model should include a constant term, False otherwise

    Returns
    -------
    None
    """
    p=order[0]
    d=order[1]
    q=order[2]

    if add_constant:
        const_coeff='c'
    else:
        const_coeff='n'
    

    if d!=0:
        x=x.diff(d)
        
    results = arimafromlib(x, order=(p,0,q), trend=const_coeff).fit()

    # fixing the wrong values in the ARIMA returned object
    results.model.order = (p,d,q)

    # fixing the wrong residuals and fittedvalues in the ARIMA returned object
    results.resid[:np.max(results.model.order)] = np.nan
    results.fittedvalues[:np.max(results.model.order)] = np.nan

    return results

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
