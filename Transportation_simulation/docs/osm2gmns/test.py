import osm2gmns as og

id = 1128379
target = 'test.osm'

og.downloadOSMData(id, target)

net = og.getNetFromFile(target)
og.show(net)