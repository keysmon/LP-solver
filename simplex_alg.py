# Hang Ruan
# V00923058
import numpy as np
import sys

non_basic = []
basic = []
c = []
def main():
    global non_basic
    global basic
    global c
    
    max_val = 0
    #zeta,row_values,coeff = file_reader('input2.txt')
    #zeta,row_values,coeff = file_reader("netlib_klein2.txt")
   
    coeff,row_values,zeta = input_reader()
    
   

    c = zeta
    n_var = len(zeta)
    n_constraint = len(row_values)
    zeta = zeta + [0]*n_constraint
    for i in range(n_var):
        non_basic.append(i)
    for i in range(n_constraint):
        basic.append(i+n_var)
        
    # Convert matrixs to numpy array
    coeff = np.array(coeff).astype(float)
    row_values = np.array(row_values).astype(float)
    zeta = np.array(zeta).astype(float)
    
    # Convert dictionary form to equational form
    identity = np.identity(n_constraint)
    simplex_table = np.hstack((coeff,identity))
   
    
    # dual simplex method to solve intially infeasible dictionary
    if not check_feasible(simplex_table,row_values,zeta):
        simplex_table,row_values,zeta,const_of_obj = dual_simplex(simplex_table,row_values,zeta)
        max_val += const_of_obj
#    print("----------------\n")
    
    cycle = 0
    use_bland_rule = 0
    while not check_optimal(simplex_table,row_values,zeta):
        #print("\n----")
        if not check_bounded(simplex_table,row_values,zeta):
            print("Unbounded")
            return
            
       
        # steepest edge rule to select enter and leaving var, bland rule if cycle
        if (use_bland_rule == 1):
            enter_index,leave_index = bland_rule(simplex_table,row_values,zeta)
            #print("Bland_rule: Enter index:",enter_index,"\tLeaving index:",leave_index)
            use_bland_rule = 0
        else:
        
            enter_index,leave_index = bland_rule(simplex_table,row_values,zeta)
            #print("Steepest_edge: Enter index:",enter_index,"\tLeaving index:",leave_index)
        if enter_index == -1:
            break
        
        
       
        # pivot operation
        simplex_table,row_values,zeta,progress = pivot(simplex_table,row_values,zeta,enter_index,leave_index)
          
        # 3 occurance of degenerate triger a blands rule
        max_val += progress
        if (progress == 0):
            if cycle == 3:
                use_bland_rule = 1
            else:
                cycle += 1
                
        
    # optimal output
    print("-----------------")
    print("optimal")
    print(round(max_val,6))
    for x in range(n_var):
        if x in basic:
            for row,row_value in zip(simplex_table,row_values):
                if  row[x] == 1:
                    print(round(row_value,6),end = " ")
        else:
            print("0",end=" ")
    print("\n-----------------\n")
    return
    
    
    
    
    
#
def dual_simplex(simplex_table,row_values,zeta):

    global basic,non_basic
    #print("\n----------------\ninfeasible to feasible")
    
    # Primal
    new_zeta = np.zeros(len(zeta))
    #primal = simplex_table[:,:len(non_basic)]
    primal = simplex_table.copy()
    new_row_values = row_values.copy()
    new_basic = basic.copy()
    new_non_basic = non_basic.copy()
    # loop while primal is not feasible
    while not check_feasible(primal,new_row_values,new_zeta):
        # identity a pair of entering and leaving var.
        unbounded_check = True
        
        #choosing enter and leaving var
        for row,row_value in zip(primal,new_row_values):
            if row_value < 0:
                # find a leave var
                for b in new_basic:
                    if row[b] == 1:
                        leave_var = b
                        # find a enter var
                        for n in new_non_basic:
                            if row[n] < 0:
                                enter_var = n
                                unbounded_check = False
                                break
                        break
                break
        if unbounded_check:
            print("infeasible")
            exit()
        #print("Enter/leave",enter_var,leave_var)
        
        #pivot operation
        primal,new_row_values,new_zeta,whatever = pivot(primal,new_row_values,new_zeta,enter_var,leave_var)
        
        #update basic and non_basic var
        new_basic.remove(leave_var)
        new_basic.append(enter_var)
        new_non_basic.remove(enter_var)
        new_non_basic.append(leave_var)
        new_basic.sort()
        new_non_basic.sort()
        
    # update new zeta
    const_of_obj = 0
    for x_index in range(len(zeta)):
        if x_index in new_basic:
            coeff = zeta[x_index]
            for row,row_value in zip(primal,new_row_values):
                if row[x_index] == 1:
                    
                    zeta -= row*coeff
                    const_of_obj += row_value*coeff
                    break
        
    basic = new_basic
    non_basic = new_non_basic
    
    return primal,new_row_values,zeta,const_of_obj


# taking all possible enter and leave var, find the pivot that best fit for the gradient of obj
def best_fit_finder(comparsion_list,simplex_table,row_values,zeta):
    global c
    
    # find X(old)
    x_old = []
    for x_index in range(len(non_basic)):
        if x_index in non_basic:
            x_old.append(0)
        else:
            x_old.append(basic.index(x_index))
    
    x_news = []
    # find X(new)s by performing all possible pivot
    for index in range(len(comparsion_list)):
        x_new = []
        enter_var = index
        leave_var = comparsion_list[index]
        
        if leave_var != -1:
            #print("  ",enter_var,leave_var)
            
            new_basic = basic.copy()
            new_non_basic = non_basic.copy()
            new_simplex_table = simplex_table.copy()
            new_row_values = row_values.copy()
            new_zeta = zeta.copy()
            new_simplex_table,new_row_values,new_zeta,whatever = pivot(new_simplex_table,new_row_values,new_zeta, enter_var,leave_var)

            new_basic.remove(leave_var)
            new_basic.append(enter_var)
            new_non_basic.remove(enter_var)
            new_non_basic.append(leave_var)
            new_basic.sort()
            new_non_basic.sort()
            
            
            for x_index in range(len(non_basic)):
                if x_index in new_non_basic:
                    x_new.append(0)
                else:
                    for row,row_value in zip(new_simplex_table,new_row_values):
                        if row[x_index] == 1:
                            x_new.append(row_value)
            x_news.append(x_new)
            
            
    # x_news has all the new x values after the pivot
    # we find the x_new maximize C.tranpose.dot(X(new)- X(old)) to find the best pivot
    growth = []
    for x_new in x_news:
        x_change = diff(x_new , x_old)
        growth.append(dot_product(c,x_change))
    order = growth.index(max(growth))
    
    for temp in range(len(comparsion_list)):
        if comparsion_list[temp] != -1:
            if order == 0:
                enter_var = temp
            else:
                order -= 1

    leave_var = comparsion_list[enter_var]
    return enter_var, leave_var
    
# find li1 - li2
def diff(li1, li2):
    diff = []
    for x,y in zip(li1,li2):
        diff.append(x-y)
    return diff

# dot product of two list
def dot_product(l1,l2):
    return sum([x*y for x,y in zip(l1,l2)])

# steepest take the pivot that max C.tranpose dot [X(new) - X(old)]
def steepest_edge_rule(simplex_table,row_values,zeta):
    bool = True
    
    # find all possible enter and leave var.
    enter_leave_list = np.zeros(len(zeta)).astype(int)
    for n in range(len(zeta)):
        if n in non_basic and zeta[n]>0:
            # find the leaving variable if a enter var is chosen
            leave_comparsion_list = np.zeros(0)
            bool2 = True

            for row,row_value in zip(simplex_table,row_values):
                if row[n]>0:
                    leave_comparsion_list = np.append(leave_comparsion_list,row_value/row[n])
                    bool2 = False
                else:
                    leave_comparsion_list=  np.append(leave_comparsion_list,float('inf'))
            
                        
                    
            leave_order = np.argmin(leave_comparsion_list)
            for basic_index in basic:
                if simplex_table[leave_order,basic_index] == 1:
                    enter_leave_list[n] = basic_index
                    bool = False
            if bool2:
                enter_leave_list[n] = -1
            
            
        else:
            enter_leave_list[n] = -1
    if bool :
        return -1,-1
    #print(enter_leave_list)
    enter_var,leave_var =  best_fit_finder(enter_leave_list,simplex_table,row_values,zeta)
    
    basic.remove(leave_var)
    basic.append(enter_var)
    non_basic.remove(enter_var)
    non_basic.append(leave_var)
    basic.sort()
    non_basic.sort()
    return enter_var,leave_var

# return the index of entering var and leaving var. Return -1,-1 if no pivotting needed
def bland_rule(simplex_table,row_values,zeta):
    global basic,non_basic
    for non_basic_index in non_basic:
        if zeta[non_basic_index] > 0: # positive coeff of non_basic_var
            comparsion_list = np.empty(0)
            
            for row,row_value in zip(simplex_table,row_values):
                if row[non_basic_index]>0:
                    comparsion_list=  np.append(comparsion_list,row_value/row[non_basic_index])
                else:
                    comparsion_list = np.append(comparsion_list,float('inf'))
            
            leave_order = np.argmin(comparsion_list)
            for basic_index in basic:
                if simplex_table[leave_order,basic_index] == 1:
                    basic.remove(basic_index)
                    basic.append(non_basic_index)
                    non_basic.remove(non_basic_index)
                    non_basic.append(basic_index)
                    basic.sort()
                    non_basic.sort()
                    return non_basic_index, basic_index
    return -1, -1


def check_optimal(simplex_table,row_values,zeta):
    for obj_coeff in zeta:
        if obj_coeff > 0:
            return False
    return True
    
    
def check_bounded(simplex_table,row_values,zeta):
    enter_index = -1
    for enter_coeff in zeta:
        enter_index +=1
        if enter_coeff > 0 :
            for row in simplex_table:
                if row[enter_index]>0:
                    return True
    return False
    
def check_feasible(simplex_table,row_values,zeta):
    for row_value in row_values:
        if row_value < 0:
            return False
    return True
    
'''
    
def file_reader(file_name):
    #get obj function
    df = pandas.read_csv(file_name,sep = '\s+',header = None,nrows = 1)
    zeta = df.iloc[0].tolist()
    
    # get the rest of data from txt
    df = pandas.read_csv(file_name,sep = '\s+',header = None,index_col = False,skiprows = 1)
    row_values = df.iloc[:,-1].tolist()
    df.drop(columns = df.columns[-1],inplace = True)
    coeff = df.values.tolist()
    return zeta,row_values,coeff
  '''
  
def input_reader():
    input = sys.stdin.read()
    input = input.split("\n")
    #get zeta
    zeta = list(map(float,input[0].split()))
    input = input[1:]
    while "" in input:
        input.remove("")
        
    coeff = []
    row_values = []
        
    for temp in input:
        temp = list(map(float,temp.split()))
        coeff.append(temp[:-1])
        row_values.append(temp[-1])
    
    return coeff,row_values,zeta
    
    
def pivot(simplex_table,row_values,zeta,enter_index,leave_index):
    global basic,non_basic
    
    # locate the row for the leaving variable
    row_index = 0
    for row in simplex_table:
        if row[leave_index] == 1:
            break
        else:
            row_index += 1
    
    row_value = row_values[row_index]
    coeff_of_leave = simplex_table[row_index,enter_index]
    
    # update the row containing the old basic
    simplex_table[row_index]/= coeff_of_leave
    row_values[row_index] /= coeff_of_leave
    
    # update simplex table
    temp = -1
    for row,row_constant in zip(simplex_table,row_values):
        temp += 1
        if (temp == row_index):
            continue
        row_values[temp] -= row[enter_index] * row_values[row_index]
        row -= row[enter_index]*simplex_table[row_index]
        
    # update zeta
    progress = zeta[enter_index] * row_values[row_index]
    zeta -= zeta[enter_index]*simplex_table[row_index]
    return simplex_table,row_values,zeta,progress
     
     
     
     
if __name__ == '__main__':
    main()




