"""Report TG-GATEs / DrugMatrix overlap with Tox21 (mechanism assays).

Reads the three *_keyed.csv files and prints the intersection counts on both
the full InChIKey and the salt/stereo-insensitive connectivity block, plus how
much Tox21 assay data the shared compounds actually carry.
"""
import csv, collections


def load(fn):
    return [r for r in csv.DictReader(open(fn)) if r['inchikey']]


def main():
    tox = list(csv.DictReader(open('tox21_keyed.csv')))
    assays = [c for c in tox[0].keys() if c not in ('inchikey', 'connectivity', 'std_smiles')]
    tox_full = {r['inchikey'] for r in tox}
    tox_conn = {}
    for r in tox:
        tox_conn.setdefault(r['connectivity'], r)

    tg = load('tggates_keyed.csv')
    dm = load('dm_keyed.csv')

    for label, src in [('TG-GATEs', tg), ('DrugMatrix', dm)]:
        full = {r['inchikey'] for r in src if r['inchikey'] in tox_full}
        conn = {r['connectivity'] for r in src if r['connectivity'] in tox_conn}
        counts = collections.Counter()
        for c in conn:
            tr = tox_conn[c]
            for a in assays:
                if tr[a] not in ('', None):
                    counts[a] += 1
        print(f"=== {label} ∩ Tox21 ===")
        print(f"  exact InChIKey : {len(full)}")
        print(f"  connectivity   : {len(conn)}")
        if conn:
            print(f"  avg assays labelled/compound: {sum(counts.values())/len(conn):.1f} / 12")

    tg_tox = {r['connectivity'] for r in tg if r['connectivity'] in tox_conn}
    dm_tox = {r['connectivity'] for r in dm if r['connectivity'] in tox_conn}
    tg_all = {r['connectivity'] for r in tg}
    dm_all = {r['connectivity'] for r in dm}
    print("\n=== Combined ===")
    print(f"  UNION ∩ Tox21                : {len(tg_tox | dm_tox)}")
    print(f"  in both transcriptome sources: {len(tg_all & dm_all)}  (cross-platform replicates)")
    print(f"  new from DrugMatrix          : {len(dm_tox - tg_tox)}")


if __name__ == '__main__':
    main()
