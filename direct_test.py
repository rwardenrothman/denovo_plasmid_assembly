import asyncio

from sqlmodel import Session

from AppContainer.app.app import transfer_annotations, combine_results
from AppContainer.app.db_model import PlasmidSeqRun, make_engine, PlasmidSeqRunList

if __name__ == '__main__':
    engine = make_engine()

    datapoint_ids = ['pGRO-K1986', 'pGRO-K1987']
    with Session(engine) as session:
        datapoints = [session.get(PlasmidSeqRun, pid) for pid in datapoint_ids]

    for datapoint in datapoints:
        print(datapoint)

        res = asyncio.run(transfer_annotations(datapoint,,)
        print(res)
        print('')

    run_list = PlasmidSeqRunList(runs=datapoints)
    res2 = asyncio.run(combine_results(run_list))
