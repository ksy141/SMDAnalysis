C6 = ['C12', 'H2', 'O2', 'HO2', 
      'C13', 'H3', 'O3', 'HO3',
      'C14', 'H4', 'O4', 'HO4',
      'C15', 'H5', 'O5', 'HO5',
      'C16', 'H6', 'O6', 'HO6',
      'C11', 'H1']
        
PO4 = ['P', 'O13',  'O14', 'O12', 'O11']

GL  = ['C1', 'HA', 'HB', 
       'C2', 'HS', 
       'C3', 'HX', 'HY',
       'O21', 'C21', 'O22',
       'O31', 'C31', 'O32']

C21  = ['C2%d' %n for n in range(2, 7)]
C21 += ['H%dR' %n for n in range(2, 7)]
C21 += ['H%dS' %n for n in range(2, 7)]

C22  = ['C2%d' %n for n in range(7, 11)]
C22 += ['H%dR' %n for n in range(7, 11)]
C22 += ['H%dS' %n for n in range(7, 11)]

C23  = ['C2%d' %n for n in range(11, 15)]
C23 += ['H%dR' %n for n in range(11, 15)]
C23 += ['H%dS' %n for n in range(11, 15)]

C24  = ['C2%d' %n for n in range(15, 19)]
C24 += ['H%dR' %n for n in range(15, 19)]
C24 += ['H%dS' %n for n in range(15, 19)]

C25  = ['C2%d' %n for n in range(19, 21)]
C25 += ['H%dR' %n for n in range(19, 21)]
C25 += ['H%dS' %n for n in range(19, 21)]
C25 += ['H20T']


C31  = ['C3%d' %n for n in range(2, 7)]
C31 += ['H%dX' %n for n in range(2, 7)]
C31 += ['H%dY' %n for n in range(2, 7)]

C32  = ['C3%d' %n for n in range(7, 11)]
C32 += ['H%dX' %n for n in range(7, 11)]
C32 += ['H%dY' %n for n in range(7, 11)]

C33  = ['C3%d' %n for n in range(11, 15)]
C33 += ['H%dX' %n for n in range(11, 15)]
C33 += ['H%dY' %n for n in range(11, 15)]

C34  = ['C3%d' %n for n in range(15, 19)]
C34 += ['H%dX' %n for n in range(15, 19)]
C34 += ['H%dY' %n for n in range(15, 19)]
C34 += ['H18Z']


group = {'C6':  C6,
         'PO4': PO4,
         'GL':  GL,
         'C21': C21,
         'C22': C22,
         'C23': C23,
         'C24': C24,
         'C31': C31,
         'C32': C32,
         'C33': C33,
         'C34': C34}

