## Finding and plotting Kangerlussuaq network - sanity check
## 1 Apr 2018 - EHU

kanger_termcoords = [(490411.71852353757, -2291155.9238763503),
 (491991.2048301102, -2289843.9028170863),
 (491201.46167682385, -2290440.2760258424),
 (489621.9753702513, -2291871.5717268577),
 (489110.9650945954, -2293004.6808234947)]

## Testing functionality of network_selection for finding and filtering
kanger_lines = {}
for j in range(len(kanger_termcoords)):
        k = kanger_termcoords[j]
        line_coords, width = Trace_wWidth(k[0], k[1], trace_up=True)
        xyw = [(line_coords[n][0], line_coords[n][1], width[n]) for n in range(len(line_coords))]
        kanger_lines[j] = (xyw)

filtered_tribs = FilterMainTributaries(kanger_lines)

## Checking that tributaries look reasonable
#plt.figure()
#plt.contourf(X, Y, S, cmap='Purples_r')
#for k in range(len(filtered_tribs)):
#    line = np.asarray(filtered_tribs[k])
#    plt.plot(line[:,0], line[:,1])
#plt.show()

## Writing a CSV file of branches traced up across terminus
## Note: WriteNetwork applies Trace_wWidth and FilterMainTributaries itself.  If you're confident in the terminus points you've chosen, no need to run those functions separately.
#WriteNetwork(kanger_termcoords[:,0], kanger_termcoords[:,1], trace_up=True, output_name='kangerlussuaq-network-w_width.csv')

#Reading in with new features of Flowline_CSV and make_full_lines
kangercoords_0, kangercoords_1, kangercoords_2, kangercoords_3, kangercoords_4 = Flowline_CSV('kangerlussuaq-network-w_width.csv', 5, has_width=True, flip_order=False)
kanger_0 = Branch(coords=kangercoords_0, index=0, order=0)
kanger_1 = Branch(coords=kangercoords_1, index=1, order=1, flows_to=0, intersect=174)
kanger_2 = Branch(coords=kangercoords_2, index=2, order=1, flows_to=0, intersect=191)
kanger_3 = Branch(coords=kangercoords_3, index=3, order=1, flows_to=0, intersect=146)
kanger_4 = Branch(coords=kangercoords_4, index=4, order=1, flows_to=0, intersect=61)
kanger_branches = (kanger_0, kanger_1, kanger_3, kanger_4)
Kanger = PlasticNetwork(name='Kangerlussuaq', init_type='Branch', branches=kanger_branches, main_terminus=kangercoords_0[0])
Kanger.make_full_lines()


plt.figure()
plt.contourf(X, Y, S, cmap='Blues_r')
for fl in Kanger.flowlines:
    plt.plot(fl.coords[:,0], fl.coords[:,1], marker='*')
plt.show()