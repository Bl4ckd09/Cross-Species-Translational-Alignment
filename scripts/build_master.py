import csv, collections
def load(fn): return [r for r in csv.DictReader(open(fn)) if r['inchikey']]

tox=list(csv.DictReader(open('tox21_keyed.csv')))
assays=[c for c in tox[0].keys() if c not in ('inchikey','connectivity','std_smiles')]
tox_conn={}
for r in tox: tox_conn.setdefault(r['connectivity'], r)   # first wins

tg=load('tggates_keyed.csv'); dm=load('dm_keyed.csv')

# index transcriptome sources by connectivity
tg_by=collections.defaultdict(list); dm_by=collections.defaultdict(list)
for r in tg: tg_by[r['connectivity']].append(r)
for r in dm: dm_by[r['connectivity']].append(r)

conns = sorted(set(tg_by)|set(dm_by))
def yn(b): return 'Y' if b else ''

fields=(['compound_name','inchikey','connectivity',
         'in_tggates','tggates_liver','tggates_kidney','in_drugmatrix',
         'cross_platform_replicate','in_tox21','tox21_assays_labelled','tox21_assays_active']
        + assays)

rows=[]
for c in conns:
    in_tg=c in tg_by; in_dm=c in dm_by
    # pick a display name: prefer TG-GATEs (drug names), else DrugMatrix
    if in_tg:
        rec=tg_by[c][0]; name=rec['name']; ik=rec['inchikey']
        liver=any(x['liver']=='True' for x in tg_by[c])
        kidney=any(x.get('kidney')=='True' for x in tg_by[c])
    else:
        rec=dm_by[c][0]; name=rec['name'].title(); ik=rec['inchikey']; liver=kidney=False
    tr=tox_conn.get(c)
    if tr:
        vals={a:tr[a] for a in assays}
        labelled=sum(1 for a in assays if tr[a] not in ('',None))
        active=sum(1 for a in assays if tr[a]=='1')
    else:
        vals={a:'' for a in assays}; labelled=active=''
    row={'compound_name':name,'inchikey':ik,'connectivity':c,
         'in_tggates':yn(in_tg),'tggates_liver':yn(in_tg and liver),
         'tggates_kidney':yn(in_tg and kidney),'in_drugmatrix':yn(in_dm),
         'cross_platform_replicate':yn(in_tg and in_dm),
         'in_tox21':yn(tr is not None),
         'tox21_assays_labelled':labelled,'tox21_assays_active':active}
    row.update(vals); rows.append(row)

rows.sort(key=lambda r:(r['in_tox21']=='', r['compound_name'].lower()))
with open('master_cohort.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

# summary
n=len(rows)
tg_n=sum(1 for r in rows if r['in_tggates'])
dm_n=sum(1 for r in rows if r['in_drugmatrix'])
tox_n=sum(1 for r in rows if r['in_tox21'])
xp=sum(1 for r in rows if r['cross_platform_replicate'])
tox_and_transc=sum(1 for r in rows if r['in_tox21'])
liver_tox=sum(1 for r in rows if r['tggates_liver'] and r['in_tox21'])
print(f"master rows (unique compounds, TG-GATEs ∪ DrugMatrix): {n}")
print(f"  in TG-GATEs: {tg_n}   in DrugMatrix: {dm_n}   cross-platform (both): {xp}")
print(f"  with Tox21 mechanism labels: {tox_n}")
print(f"  TG-GATEs liver + Tox21: {liver_tox}")
