
GL  = ['C1',  'HA',  'HB',
       'O11', 'C11', 'O12',
       'C2',  'HS', 
       'O21', 'C21', 'O22', 
       'C3',  'HX',  'HY',
       'O31', 'C31', 'O32']

C11  = ['C1%d' %n for n in range(2, 7)]
C11 += ['H%dA' %n for n in range(2, 7)]
C11 += ['H%dB' %n for n in range(2, 7)]

C12  = ['C1%d' %n for n in range(7, 11)]
C12 += ['H%dA' %n for n in range(7, 11)]
C12 += ['H%dB' %n for n in range(7, 11)]

C13  = ['C1%d' %n for n in range(11, 15)]
C13 += ['H%dA' %n for n in range(11, 15)]
C13 += ['H%dB' %n for n in range(11, 15)]

C14  = ['C1%d' %n for n in range(15, 19)]
C14 += ['H%dA' %n for n in range(15, 19)]
C14 += ['H%dB' %n for n in range(15, 19)]
C14 += ['H18C']

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
C24 += ['H18T']

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


group = {'GL':  GL,
         'C11': C11,
         'C12': C12,
         'C13': C13,
         'C14': C14,
         'C21': C21,
         'C22': C22,
         'C23': C23,
         'C24': C24,
         'C31': C31,
         'C32': C32,
         'C33': C33,
         'C34': C34}

