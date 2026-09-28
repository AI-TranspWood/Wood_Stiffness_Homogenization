"""Tkinter front end for hom_tw.py: pick a wood, edit inputs, homogenize.

Every numeric field accepts one value, a list '500, 600, 700' or a MATLAB-style
range '500:50:800'; all combinations are computed for native, infiltrated and
delignified + infiltrated wood. The field keys are the input-file keys, so
'Save input' / 'Load input' write and read the same JSON as `python hom_tw.py`.
"""
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import hom_tw as h

root = tk.Tk()
root.title('Transparent wood - stiffness homogenization')
V = {}                                              # input key -> tk variable
res = dict(rows=[], varied=[])                      # last run_batch result
COLOR = {1: '#8c5a2b', 2: '#1f77b4', 3: '#2ca02c'}  # native / infiltrated / delignified + infiltrated
BARS = ['E_R', 'E_T', 'E_L', 'G_TL', 'G_RL', 'G_RT', 'M_cw']
WOOD_KEYS = ('density', 'CI', 'MFA', 'cellulose', 'hemicellulose', 'lignin', 'extractives',
             'vessel_pct', 'ray_pct', 'L_EW', 'W2_EW', 'L_LW', 'W2_LW')


def field(parent, row, key, width=12):
    V[key] = tk.StringVar()
    ttk.Label(parent, text=h.INPUTS.get(key, key)).grid(row=row, column=0, sticky='w', padx=4, pady=1)
    ttk.Entry(parent, textvariable=V[key], width=width).grid(row=row, column=1, sticky='w', padx=4)


def fmt(x):
    return 'NaN' if not np.isfinite(x) else f'{x:g}'


def show(key, v):                                   # value(s) from defaults / input file -> field text
    if key in ('wood', 'polymer', 'cellulose_material'):
        V[key].set(v)
    elif key == 'swelling':
        V[key].set(int(bool(v)))
    elif key in V:
        V[key].set(v if isinstance(v, str) else ', '.join(fmt(x) for x in h.values(v)))


def fill_from_wood(*_):
    f = h.flat_defaults(V['wood'].get())
    for k in WOOD_KEYS:
        show(k, f[k])
    wood_info()


def wood_info():
    w = h.WOODS[V['wood'].get()]
    info.set(f"{w.get('latin', '')}, {w.get('type', '')}; density {w['density_kg_m3']['min']}-"
             f"{w['density_kg_m3']['max']} kg/m³, MFA {w['MFA_deg']['min']}-{w['MFA_deg']['max']}°"
             + (f"\nEW/LW geometry: {w['EWLW_surrogate']} surrogate" if 'EWLW_surrogate' in w else ''))


def fill_polymer(*_):
    if V['polymer'].get() in h.POLYMERS:
        for k, v in zip(h.POLY, h.POLYMERS[V['polymer'].get()]):
            V[k].set(fmt(v))


def inputs():
    """Current fields as an input dict (numbers where single, text where list/range)."""
    d = dict(wood=V['wood'].get(), states=[1, 2, 3], polymer=V['polymer'].get(),
             cellulose_material=V['cellulose_material'].get(), swelling=bool(V['swelling'].get()))
    for k in list(h.INPUTS) + ['MFA']:
        t = V[k].get().strip()
        try:
            d[k] = float(t) if np.isfinite(float(t)) else None     # NaN -> JSON null
        except ValueError:
            d[k] = t
    return d


def save_input():
    fn = filedialog.asksaveasfilename(defaultextension='.json', filetypes=[('input file', '*.json')])
    if fn:
        with open(fn, 'w', encoding='utf8') as fh:
            json.dump(inputs(), fh, indent=1, ensure_ascii=False)
        status.set(f'saved {fn}')


def load_input(fn=None):
    fn = fn or filedialog.askopenfilename(filetypes=[('input file', '*.json')])
    if not fn:
        return
    with open(fn, encoding='utf8') as fh:
        inp = json.load(fh)
    V['wood'].set(inp.get('wood', V['wood'].get()))
    f = h.flat_defaults(V['wood'].get())
    if inp.get('polymer') in h.POLYMERS and 'E_poly' not in inp:
        f['E_poly'], f['nu_poly'] = h.POLYMERS[inp['polymer']]
    f.update(inp)
    for k, v in f.items():
        show(k, v)
    wood_info()
    status.set(f'loaded {fn}')


def run():
    try:
        b = h.run_batch(inputs(), lambda n, N: (status.set(f'run {n}/{N}'), root.update_idletasks()))
    except ValueError as e:
        messagebox.showerror('Input', str(e))
        return
    status.set(f"{b['n']} run(s), {len(b['rows'])} result row(s)"
               + (f", {len(b['errors'])} failed" if b['errors'] else ''))
    if b['errors']:
        messagebox.showwarning('Some runs failed', '\n'.join(b['errors'][:10]))
    res.clear()
    res.update(b)
    show_table()
    xsel['values'] = [label(k) for k in res['varied']]
    xsel.set(xsel['values'][0] if res['varied'] else '')
    build_selectors()
    if res['rows']:
        tree.selection_set(tree.get_children()[0])


def label(k):
    return {'state': 'state', 'MFA': 'MFA [°]'}.get(k) or h.INPUTS[k]


def show_table():
    ins = ['state'] + [k for k in res['varied'] if k != 'MFA'] + ['MFA']
    tree.delete(*tree.get_children())
    tree['columns'] = ins + list(h.OUTPUTS)
    for c in tree['columns']:
        tree.heading(c, text=label(c) if c in ins else h.OUTPUTS[c])
        tree.column(c, width={'state': 150}.get(c, 95 if c in ins else 72), anchor='e', stretch=False)
    for j, (d, r, i) in enumerate(res['rows']):
        tree.insert('', 'end', iid=str(j), tags=(str(d['state']),),
                    values=[h.STATES[d['state']]] + [fmt(d[k]) for k in ins[1:]]
                    + [f'{r[k][i]:.4g}' for k in h.OUTPUTS])
    for st, c in COLOR.items():
        tree.tag_configure(str(st), foreground=c)


def show_row(*_):
    sel = tree.selection()
    if not sel:
        return
    d, r, i = res['rows'][int(sel[0])]
    v, s = r['vol'], r['vol']['swelling']
    txt.delete('1.0', 'end')
    txt.insert('end', h.STATES[d['state']] + ''.join(f', {label(k)} = {fmt(d[k])}' for k in res['varied']) + '\n'
               + f"f_bcw (clearwood) = {v['fclw']['bcw']:.4f}   EW/LW = {v['EW']:.3f}/{v['LW']:.3f}   "
               + 'cell wall: ' + ', '.join(f'{k} {x:.3f}' for k, x in v['fcw'].items()) + '\n')
    if s['pct']:
        txt.insert('end', f"swelling {s['pct']:g} %: wall volume x{s['k_wall']:.4f}, fibre lumen porosity "
                          f"{s['porosity_fib']:.4f} -> {s['porosity_fib_swollen']:.4f}\n")
    for name, C in (('clearwood', r['C_clearwood'][i]), ('cell wall', r['C_cellwall'][i])):
        txt.insert('end', f'\nC {name}, MFA = {fmt(d["MFA"])}° [GPa, Kelvin-Mandel]\n')
        txt.insert('end', '\n'.join(' '.join(f'{x:9.3f}' for x in row) for row in C) + '\n')


SLICE = {}                                          # other varied input -> combobox of its values


def xkey():
    return next((k for k in res['varied'] if label(k) == xsel.get()), res['varied'][0])


def build_selectors(*_):
    """No batch: bars. One varied input: y selector. Several: y, x and values of the others."""
    for w in slf.winfo_children():
        w.destroy()
    SLICE.clear()
    nb = len(res['varied'])
    for w, on in ((ylab, nb > 0), (ysel, nb > 0), (xlab, nb > 1), (xsel, nb > 1)):
        w.grid() if on else w.grid_remove()
    if nb > 1:
        for j, k in enumerate(k for k in res['varied'] if k != xkey()):
            vals = sorted({d[k] for d, _, _ in res['rows']})
            ttk.Label(slf, text=label(k) + ' =').grid(row=0, column=2 * j, padx=(8, 2))
            SLICE[k] = ttk.Combobox(slf, values=[fmt(x) for x in vals], state='readonly', width=8)
            SLICE[k].set(fmt(vals[0]))
            SLICE[k].grid(row=0, column=2 * j + 1)
            SLICE[k].bind('<<ComboboxSelected>>', plot)
    plot()


def plot(*_):
    ax.clear()
    rows = res['rows']
    if not rows:
        canvas.draw()
        return
    if not res['varied']:                           # single set: grouped bars per state
        w = 0.8 / 3
        for d, r, i in rows:
            ax.bar(np.arange(len(BARS)) + (d['state'] - 2) * w, [r[k][i] for k in BARS], w,
                   color=COLOR[d['state']], label=h.STATES[d['state']])
        ax.set_xticks(range(len(BARS)), BARS)
        ax.set_ylabel('[GPa]')
        ax.grid(alpha=0.3, axis='y')
    else:                                           # x-y plot, one line per state
        key = next(k for k, lab in h.OUTPUTS.items() if lab == ysel.get())
        x = xkey()
        fixed = {k: float(cb.get()) for k, cb in SLICE.items()}
        xs = sorted({d[x] for d, _, _ in rows})
        for st in h.STATES:
            pts = sorted((d[x], r[key][i]) for d, r, i in rows if d['state'] == st and all(
                np.isclose(d[k], v) for k, v in fixed.items() if not (st == 1 and k in h.POLY)))
            if st == 1 and x in h.POLY and pts:     # native does not depend on the polymer
                ax.plot([xs[0], xs[-1]], [pts[0][1]] * 2, '--', color=COLOR[st], label=h.STATES[st])
            elif pts:
                ax.plot(*zip(*pts), 'o-', ms=3, color=COLOR[st], label=h.STATES[st])
        ax.set_xlabel(label(x))
        ax.set_ylabel(ysel.get())
        ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    canvas.draw()


def save_csv():
    if res['rows']:
        fn = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[('CSV', '*.csv')])
        if fn:
            h.write_csv(res, fn)
            status.set(f'saved {fn}')


# ------------------------------------------------------------------ layout
left = ttk.Frame(root, padding=6)
left.grid(row=0, column=0, sticky='ns')
right = ttk.Frame(root, padding=6)
right.grid(row=0, column=1, sticky='nsew')
root.columnconfigure(1, weight=1)
root.rowconfigure(0, weight=1)

V['wood'] = tk.StringVar(value='Birch' if 'Birch' in h.WOODS else next(iter(h.WOODS)))
ttk.Label(left, text='Wood').grid(row=0, column=0, sticky='w', padx=4)
cb = ttk.Combobox(left, textvariable=V['wood'], values=list(h.WOODS), state='readonly', width=14)
cb.grid(row=0, column=1, sticky='w', padx=4)
cb.bind('<<ComboboxSelected>>', fill_from_wood)
info = tk.StringVar()
ttk.Label(left, textvariable=info, foreground='gray', wraplength=300).grid(row=1, column=0, columnspan=2, sticky='w')
ttk.Label(left, text="numeric fields: '600', '550, 650' or range '500:50:800'\n"
                     "-> all combinations are computed", foreground='#1a5fb4').grid(
    row=2, column=0, columnspan=2, sticky='w', pady=(4, 2))

ba = ttk.Frame(left)
ba.grid(row=3, column=0, columnspan=2, sticky='ew')
field(ba, 0, 'density')
field(ba, 1, 'moisture')
ttk.Label(left, text='native, infiltrated and delignified + infiltrated\nare always computed (brown / blue / green)',
          foreground='gray').grid(row=4, column=0, columnspan=2, sticky='w', pady=4)

po = ttk.LabelFrame(left, text='Polymer and interface', padding=4)
po.grid(row=5, column=0, columnspan=2, sticky='ew', pady=4)
V['polymer'] = tk.StringVar(value='HEMA')
ttk.Label(po, text='Polymer').grid(row=0, column=0, sticky='w', padx=4)
pc = ttk.Combobox(po, textvariable=V['polymer'], values=list(h.POLYMERS) + ['custom'], state='readonly', width=11)
pc.grid(row=0, column=1, sticky='w', padx=4)
pc.bind('<<ComboboxSelected>>', fill_polymer)
for i, k in enumerate(('E_poly', 'nu_poly', 'IF_alpha', 'IF_beta'), 1):
    field(po, i, k)

sw = ttk.LabelFrame(left, text='Swelling', padding=4)
sw.grid(row=6, column=0, columnspan=2, sticky='ew', pady=4)
V['swelling'] = tk.IntVar(value=0)
ttk.Checkbutton(sw, text='cell walls swell', variable=V['swelling']).grid(row=0, column=0, columnspan=2, sticky='w')
field(sw, 1, 'swelling_pct')
ttk.Label(sw, text='ML+S2 wall thickens into the lumen; the extra\nwall volume is pore space (polymer when infiltrated)',
          foreground='gray').grid(row=2, column=0, columnspan=2, sticky='w')

adv = ttk.LabelFrame(left, text='Advanced', padding=4)


def toggle_adv():
    on = not adv.winfo_ismapped()
    adv.grid(row=8, column=0, columnspan=2, sticky='ew') if on else adv.grid_remove()
    adv_btn.config(text=('▾' if on else '▸') + ' Advanced')


adv_btn = ttk.Button(left, text='▸ Advanced', command=toggle_adv)
adv_btn.grid(row=7, column=0, columnspan=2, sticky='w', pady=2)
V['MFA'] = tk.StringVar()
ttk.Label(adv, text='MFA [°] (list per run)').grid(row=0, column=0, sticky='w', padx=4)
ttk.Entry(adv, textvariable=V['MFA'], width=16).grid(row=0, column=1, sticky='w', padx=4)
for i, k in enumerate(('CI', 'cellulose', 'hemicellulose', 'lignin', 'extractives', 'vessel_pct', 'ray_pct',
                       'L_EW', 'W2_EW', 'L_LW', 'W2_LW', 'cell_aspect_ratio', 'superellipse_n',
                       'n_families', 'microfibril_slenderness', 'tolerance'), 1):
    field(adv, i, k)
V['cellulose_material'] = tk.StringVar()
ttk.Label(adv, text='Crystalline cellulose').grid(row=20, column=0, sticky='w', padx=4)
ttk.Combobox(adv, textvariable=V['cellulose_material'], values=list(h.CRYCEL), state='readonly',
             width=11).grid(row=20, column=1, sticky='w', padx=4)

bf = ttk.Frame(left)
bf.grid(row=9, column=0, columnspan=2, pady=8)
for j, (t, cmd) in enumerate((('Run', run), ('Save CSV', save_csv), ('Save input…', save_input),
                              ('Load input…', load_input))):
    ttk.Button(bf, text=t, command=cmd).grid(row=j // 2, column=j % 2, padx=4, pady=2, sticky='ew')
status = tk.StringVar()
ttk.Label(left, textvariable=status, foreground='gray', wraplength=300).grid(row=10, column=0, columnspan=2, sticky='w')

tf = ttk.Frame(right)
tf.grid(row=0, column=0, sticky='nsew')
tree = ttk.Treeview(tf, show='headings', height=10, selectmode='browse')
ys = ttk.Scrollbar(tf, orient='vertical', command=tree.yview)
xs_ = ttk.Scrollbar(tf, orient='horizontal', command=tree.xview)
tree.configure(yscrollcommand=ys.set, xscrollcommand=xs_.set)
tree.grid(row=0, column=0, sticky='nsew')
ys.grid(row=0, column=1, sticky='ns')
xs_.grid(row=1, column=0, sticky='ew')
tf.columnconfigure(0, weight=1)
tree.bind('<<TreeviewSelect>>', show_row)

pf = ttk.Frame(right)
pf.grid(row=1, column=0, sticky='w', pady=(6, 0))
ylab = ttk.Label(pf, text='plot')
ylab.grid(row=0, column=0)
ysel = ttk.Combobox(pf, values=list(h.OUTPUTS.values()), state='readonly', width=16)
ysel.set('E_L [GPa]')
ysel.grid(row=0, column=1, padx=4)
xlab = ttk.Label(pf, text='vs')
xlab.grid(row=0, column=2)
xsel = ttk.Combobox(pf, state='readonly', width=24)
xsel.grid(row=0, column=3, padx=4)
slf = ttk.Frame(pf)                                 # values of the other varied inputs
slf.grid(row=0, column=4)
ysel.bind('<<ComboboxSelected>>', plot)
xsel.bind('<<ComboboxSelected>>', build_selectors)

fig = Figure(figsize=(6, 3.2), dpi=100)
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.16)
canvas = FigureCanvasTkAgg(fig, master=right)
canvas.get_tk_widget().grid(row=2, column=0, sticky='nsew', pady=4)
txt = tk.Text(right, height=14, width=80, font=('Consolas', 9))
txt.grid(row=3, column=0, sticky='nsew')
right.columnconfigure(0, weight=1)
right.rowconfigure(2, weight=1)
right.rowconfigure(3, weight=1)

for k, v in h.flat_defaults(V['wood'].get()).items():
    show(k, v)
fill_from_wood()

if __name__ == '__main__':
    root.mainloop()
