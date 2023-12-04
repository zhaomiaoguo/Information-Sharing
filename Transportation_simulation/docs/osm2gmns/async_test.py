import osm2gmns as og
import asyncio

target = 'test.osm'
id = 1128379

async def dl():
    og.downloadOSMData(id, target)

async def main():
    await dl()
    net = og.getNetFromFile(target)
    og.show(net)

if __name__ == '__main__':
    asyncio.run(main())