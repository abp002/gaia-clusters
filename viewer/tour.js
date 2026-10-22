// The guided tour: a few stops per region. A stop names a cluster (Hunt & Reffert's name) or none for the
// whole region, and says why it is worth a look. The numbers come from the region's data (`c` is the
// cluster's entry in meta.json, `m` the region's meta, `f` formatters), so the text cannot drift from them.

export const TOURS = {
  carina: [
    { text: (c, m, f) => `Carina, Galactic longitude 285–295°: the 10° window with the most Hunt & Reffert clusters.
      ${m.recovered} of ${m.recoverable} come back, with the selection method fixed before the run.` },
    { cluster: 'IC_2602', text: (c, m, f) => `IC 2602, the Southern Pleiades, ${f.kpc(c.dist)} away: the nearest cluster here.
      So close, its ${f.num(c.n)} members spread over several degrees of sky, and HDBSCAN still finds ${f.pct(c.recovery)} of them.` },
    { cluster: 'NGC_3532', text: (c, m, f) => `NGC 3532, the Wishing Well Cluster: the richest here, with ${f.num(c.n)} members in our sample.
      Recovery ${f.pct(c.recovery)}, purity ${f.pct(c.purity)}.` },
    { cluster: 'Trumpler_16', text: (c, m, f) => `The Carina Nebula, about ${f.kpc(c.dist)} away, where stars are still forming.
      Trumpler 16 and Trumpler 14 are a few million years old and both come back; Collinder 228, beside them,
      is swallowed by a bigger group.` },
    { cluster: 'UPK_562', text: (c, m, f) => `UPK 562 misses by a hair. The group holds ${f.pct(c.recovery)} of its members,
      but they are only ${f.pct(c.purity)} of the group, and the rule asks for more than half.
      ${f.purityMisses(m)} of Carina's ${m.recoverable - m.recovered} misses fail like this, on purity.` },
    { cluster: 'Trumpler_19', text: (c, m, f) => `Trumpler 19, ${f.age(c)}: one of the oldest clusters here.
      Old clusters thin out as their stars drift away, but this one keeps ${f.num(c.n)} members, and all of them come back.` },
    { cluster: 'HSC_2393', text: (c, m, f) => `HSC 2393, ${f.kpc(c.dist)} away, the farthest cluster in the region.
      Only ${f.num(c.n)} of its stars pass our cuts out here, and it is missed: few stars in our sample are this far.` },
  ],
  auriga: [
    { text: (c, m, f) => `Auriga, Galactic longitude 170–180°: the direction away from the Galactic centre, where the field is thinnest.
      ${m.recovered} of ${m.recoverable} come back with leaf selection, chosen after the protocol's own method failed.` },
    { cluster: 'NGC_2099', text: (c, m, f) => `M37, ${f.kpc(c.dist)} away, the richest of Auriga's three Messier clusters,
      with ${f.num(c.n)} members in our sample. Recovery ${f.pct(c.recovery)}, purity ${f.pct(c.purity)}.` },
    { cluster: 'NGC_1960', text: (c, m, f) => `M36, ${f.age(c)}, much younger than M37 and M38 beside it.
      Its group takes in more field stars: purity ${f.pct(c.purity)}.` },
    { cluster: 'Stock_8', text: (c, m, f) => `Stock 8, ${f.age(c)}, one of the youngest here.
      Its group is almost all members (purity ${f.pct(c.purity)}), but it holds ${f.pct(c.recovery)} of them.` },
    { cluster: 'UBC_1281', text: (c, m, f) => `UBC 1281 is missed on the line: its members make up exactly half of the group,
      and the rule asks for more than half.` },
    { cluster: 'King_8', text: (c, m, f) => `King 8, ${f.kpc(c.dist)} away, the farthest here. ${f.num(c.n)} members pass our cuts,
      the fewest a cluster can have and still count, and it is missed.` },
  ],
};
