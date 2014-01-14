#!/usr/bin/python2.7

ufile = open("fisher_weights_normalized.txt")

text = ""
variables = []
values = []
for line in ufile.readlines():
    if "Expression" in line:
        on = False
        var = ""
        for i in range(11, len(line)):
            if on and line[i] == '"':
                on = False
            if on:
                var += line[i]
            if 'Expression=' in line[i-12:i]:
                on = True
        variables.append(var)
        
    if "Value" in line:
        on = False
        var = ""
        for i in range(6, len(line)):
            if on and line[i] == '"':
                on = False
            if on:
                var += line[i]
            if 'Value=' in line[i-6:i]:
                on = True
        values.append(abs(float(var)))
        
fisher_list = zip(values, variables)
fisher_list.sort()
fisher_list = reversed(fisher_list)

for val, var in fisher_list:
    print var, val