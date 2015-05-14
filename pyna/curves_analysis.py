import pandas as pd
import numpy as np
import re
import math
import matplotlib.pyplot as plt

def floatify(val):
    try:
        return float(val)
    except ValueError:
        return float("nan")
    
def make_key(line):
    val2 = re.sub("\)", "", line[0])
    return float(val2)

class StrandStatistics(object):
    def __init__(self, df, locations=None, times=None):
        # regularize locations to be a list
        self.initial_df = df
        if locations == None:
            locations = df.columns.values
        self.locations = locations
        if times == None:
            times = df.index.values
        self.times = times
        
        self.df = pd.DataFrame(df.loc[times, locations])
        np_arr = self.df.values.flatten()
        self.np = np_arr[~np.isnan(np_arr)]

    def _stat_obj(self, per_location, basepair_label):
        if per_location == True:
            return self.df
        else:
            return self.np
    
    def summary(self, per_location=False, basepair_label=False):
        wrap = lambda x : str(x.values) if hasattr(x, 'values') else str(x)
        summ = ""
        summ += "count: " + wrap(self.count(per_location)) + "\n"
        summ += "mean: " + wrap(self.mean(per_location)) + "\n"
        summ += "std:  " + wrap(self.std(per_location)) + "\n"
        summ += "min:  " + wrap(self.min(per_location)) + "\n"
        summ += "max:  " + wrap(self.max(per_location)) + "\n"
        return summ

    def count(self, per_location=False, basepair_label=False):
        """Count of the number of non-NaN results."""
        # TODO: this isn't quite right
        if per_location == False:
            return len(self.np)
        else:
            return self.df.count()
        
    def mean(self, per_location=False, basepair_label=False):
        stat_obj = self._stat_obj(per_location, basepair_label)
        if len(stat_obj) > 0:
            return stat_obj.mean()
        else:
            return float('nan')
    
    def std(self, per_location=False, basepair_label=False):
        stat_obj = self._stat_obj(per_location, basepair_label)
        if len(stat_obj) > 0:
            return stat_obj.std()
        else:
            return float('nan')
    
    def median(self, per_location=False, basepair_label=False):
        stat_obj = self._stat_obj(per_location, basepair_label)
        if len(stat_obj) > 0:
            return stat_obj.median()
        else:
            return float('nan')
    
    def min(self, per_location=False, basepair_label=False):
        stat_obj = self._stat_obj(per_location, basepair_label)
        return stat_obj.min()
    
    def max(self, per_location=False, basepair_label=False):
        stat_obj = self._stat_obj(per_location, basepair_label)
        return stat_obj.max()

    def hist(self, per_location=False, basepair_label=False, **kwargs):
        # TODO: per_location requires a bit more work here
        if per_location == False:
            if len(self.calc_np) > 0:
                return np.histogram(self.np, **kwargs)
            else:
                return [[],[]] # empty histogram for np return style
        else:
            hists = [np.histogram(self.df[col], **kwargs)
                     for col in self.df.columns]
            # get bins
            for hist in hists:
                if len(hist[1]) > 0:
                    bins = hist[1]
                    break
            # prepare a dictionary of histograms 
            hist_dict = {}
            for i in range(len(hists)):
                if len(hists[i][0]) == 0:
                    myhist = [0]*(len(bins)-1)
                else:
                    myhist = hists[i][0]
                hist_dict[self.df.columns.values[i]] = myhist
            # convert the dictionary to a DataFrame
            hist_df = pd.DataFrame(hist_dict).set_index(bins[:-1])
            return hist_df
    
    def __str__(self):
        return str(self.df)


def curves_style(panel):
    return panel.transpose('major_axis', 'minor_axis', 'items')

class CurvesAnalysis(object):
    # TODO: write docs ... basic idea, though is that we create a dict of
    # several pd.Panels.  One panel per "group" of data (curves makes 5 such
    # groups) One DataFrame per "measurement" (buckle, opening, minor groove
    # width) One Column per location (usually by base pair number)
    # One row per time slice
    def __init__(self, fname=None):
        self.fname = fname
        
        self.setup = {}
        self.setup['groupA'] = {
            4 : 'xdisp',
            5 : 'ydisp',
            6 : 'inclin',
            7 : 'tip',
            8 : 'ax_bend'
        }
        
        self.setup['groupB'] = {
            4 : 'shear',
            5 : 'stretch',
            6 : 'stagger',
            7 : 'buckle',
            8 : 'propel',
            9 : 'opening'
        }
        
        self.setup['groupC'] = {
            4 : 'shift',
            5 : 'slide',
            6 : 'rise',
            7 : 'tilt',
            8 : 'roll',
            9 : 'twist',
            10 : 'h-ris',
            11 : 'h-twi'
        }
        
        self.setup['groupE'] = {
            3 : "w12",
            4 : "d12",
            5 : "w21",
            6 : "d21"
        }
        
        self.group_labels = ['groupA', 'groupB', 'groupC', 'groupE']

        self.panels = {}

        self.strands = {}
        
        self.co_keys = {}
        self.prep_data = {}
        for label in self.group_labels:
            self.co_keys[label] = {}
            self.prep_data[label] = {}
            
        if self.fname is not None:
            self.read_curves_file(fname)
            
    def is_strand(self, line):
        if re.search("\(3'-5'\)", line):
            return True
        elif re.search("\(5'-3'\)", line):
            return True
        else:
            return False

    def get_strand(self, line):
        if re.search("\(3'-5'\)", line):
            key = "35"
        elif re.search("\(5'-3'\)", line):
            key = "53"
        value = re.split(": ", line)[1].rstrip()
        return (key, value)

    def is_data(self, line):
        splitter = line.split()
        try:
            is_data = re.search("[0-9\.]+\)*", splitter[0])
        except IndexError:
            is_data = False
        return is_data
        
    
    def read_curves_file(self, fname=None):
        if fname == None:
            fname = self.fname
        if fname == None:
            raise RuntimeError("No file defined for analysis!")
        f = open(fname, "r")
        category = ""
        countA = 0
        countB = 0
        countC = 0
        countD = 0
        countE = 0
        
        for line in f:
            if re.search("\(A\)", line):
                category = "A"
                countA += 1
            elif re.search("\(B\)", line):
                category = "B"
                countB += 1
            elif re.search("\(C\)", line):
                category = "C"
                countC += 1
            elif re.search("\(D\)", line):
                category = "D"
                countD += 1
            elif re.search("\(E\)", line):
                category = "E"
                countE += 1            

            if category == "A" and self.is_data(line):
                group_label = 'groupA'
                co_keys = range(1,4)
                float_data = range(4,9)
                str_data = []
                splitter = self.line_prep(line, co_keys, float_data, str_data)
                self.add_data(splitter, group_label, co_keys, float_data, str_data)
            elif category == "B" and self.is_data(line):
                group_label = 'groupB'
                co_keys = range(1,4)
                float_data = range(4,10)
                str_data = []
                splitter = self.line_prep(line, co_keys, float_data, str_data)
                self.add_data(splitter, group_label, co_keys, float_data, str_data)
            elif category == "C" and self.is_data(line):
                group_label = 'groupC'
                co_keys = range(1,4)
                float_data = range(4, 12)
                str_data = []
                splitter = self.line_prep(line, co_keys, float_data, str_data)
                self.add_data(splitter, group_label, co_keys, float_data, str_data)                
            elif category == "E" and self.is_data(line):
                group_label = 'groupE'
                co_keys = [1, 2]
                float_data = range(3, 7)
                str_data = []
                splitter = self.line_prep(line, co_keys, float_data, str_data)
                key = make_key(splitter)
                if key - math.floor(key) > 0.01: # we have an x.5 value
                    splitter.insert(1, '---')
                    splitter.insert(1, '---')
                self.add_data(splitter, group_label, co_keys, float_data, str_data)
            elif self.is_strand(line):
                (k, v) = self.get_strand(line)
                self.strands[k] = v
    
        
        for label in self.group_labels:
            dfs = {}
            for (key, name) in zip(self.setup[label].keys(), self.setup[label].values()):
                dfs[name] = pd.DataFrame(self.prep_data[label][key])
            self.panels[label] = pd.Panel(dfs)
            
        return (countA, countB, countC, countD, countE)
    
    def line_prep(self, line, co_keys, float_data, str_data):
        splitter = line.split()
        global_max = -1
        for indices in [co_keys, float_data, str_data]:
            if len(indices) > 0:
                if max(indices) > global_max:
                    global_max = max(indices)
            
        while len(splitter) < global_max+1:
            splitter.append('---')
            
        return splitter

                
    def add_data(self, splitter, group_label, co_keys, float_data, str_data):
        setup = self.setup[group_label]
        try:
            prep_data = self.prep_data[group_label]
        except KeyError:
            self.prep_data[group_label] = { }
            prep_data = self.prep_data[group_label]
        
        key = make_key(splitter)
        self.co_keys[group_label][key] = [splitter[k] for k in co_keys]
        for col in setup.keys():
            if col in float_data:
                val = floatify(splitter[col])
            elif col in str_data:
                val = str(splitter[col])
            try:
                prep_data[col][key].append(val)
            except KeyError:
                try:
                    prep_data[col][key] = [val]
                except KeyError:
                    prep_data[col] = {}
                    prep_data[col][key] = [val]

