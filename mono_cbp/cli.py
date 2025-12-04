"""
Command-line interface for mono-cbp.
"""

import argparse
import sys


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='mono-cbp: Circumbinary Planet Detection Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Run complete pipeline
    run_parser = subparsers.add_parser('run', help='Run complete pipeline')
    run_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV with eclipse and orbital parameters')
    run_parser.add_argument('--data-dir', required=True, help='Data directory')
    run_parser.add_argument('--sector-times', help='Path to sector times CSV')
    run_parser.add_argument('--output-dir', default='results', help='Output directory')
    run_parser.add_argument('--config', help='Path to configuration JSON file')
    run_parser.add_argument('--tebc', action='store_true', help='Use TEBC catalogue format (with *_2g and *_pf columns)')

    # Eclipse masking
    mask_parser = subparsers.add_parser('mask-eclipses', help='Mask eclipses only')
    mask_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV')
    mask_parser.add_argument('--data-dir', required=True, help='Data directory (files modified in-place)')
    mask_parser.add_argument('--tebc', action='store_true', help='Use TEBC catalogue format (with *_2g and *_pf columns)')

    # Transit finding
    find_parser = subparsers.add_parser('find-transits', help='Find transits only')
    find_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV with eclipse and orbital parameters')
    find_parser.add_argument('--data-dir', required=True, help='Data directory')
    find_parser.add_argument('--sector-times', help='Path to sector times CSV')
    find_parser.add_argument('--output', default='transit_events.txt', help='Output file')
    find_parser.add_argument('--threshold', type=float, help='MAD threshold')
    find_parser.add_argument('--method', choices=['cb', 'cp'], help='Detrending method')
    find_parser.add_argument('--tebc', action='store_true', help='Use TEBC catalogue format (with *_2g and *_pf columns)')

    # Model comparison
    compare_parser = subparsers.add_parser('compare-models', help='Compare models for vetting')
    compare_parser.add_argument('--event-dir', required=True, help='Event snippets directory')
    compare_parser.add_argument('--output', default='classifications.csv', help='Output file')

    # Injection-retrieval
    inject_parser = subparsers.add_parser('inject-retrieve', help='Run injection-retrieval')
    inject_parser.add_argument('--models', required=True, help='Path to transit models .npz')
    inject_parser.add_argument('--data-dir', required=True, help='Data directory')
    inject_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV with eclipse and orbital parameters')
    inject_parser.add_argument('--output', default='injection_results.csv', help='Output file')
    inject_parser.add_argument('--n-injections', type=int, default=100, help='Number of injections')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Import here to avoid slow startup
    if args.command == 'run':
        from mono_cbp import MonoCBPPipeline
        import json

        config = {}
        if hasattr(args, 'config') and args.config:
            with open(args.config) as f:
                config = json.load(f)

        pipeline = MonoCBPPipeline(
            catalogue_path=args.catalogue,
            data_dir=args.data_dir,
            output_dir=args.output_dir if hasattr(args, 'output_dir') else './results',
            sector_times_path=args.sector_times if hasattr(args, 'sector_times') else None,
            TEBC=args.tebc if hasattr(args, 'tebc') else False,
            config=config
        )

        results = pipeline.run()

        # Count detected events and candidates
        n_events = len(results.get('transit_finding', [])) if results.get('transit_finding') is not None else 0

        # Count high-confidence candidates if vetting was run
        n_candidates = 0
        if 'vetting' in results and results['vetting'] is not None:
            vetting_df = results['vetting']
            n_candidates = len(vetting_df[vetting_df['best_fit'].isin(['T', 'AT'])]) if len(vetting_df) > 0 else 0

        print(f"\nPipeline complete!")
        print(f"  Transit events detected: {n_events}")
        if 'vetting' in results:
            print(f"  High-confidence candidates: {n_candidates}")

    elif args.command == 'mask-eclipses':
        from mono_cbp import EclipseMasker
        from mono_cbp.utils import load_catalogue

        catalogue = load_catalogue(args.catalogue, TEBC=args.tebc if hasattr(args, 'tebc') else False)
        masker = EclipseMasker(catalogue, data_dir=args.data_dir)
        masker.mask_all()
        print("Eclipse masking complete!")

    elif args.command == 'find-transits':
        from mono_cbp import TransitFinder
        from mono_cbp.utils import load_catalogue

        config = {}
        if args.threshold:
            config['transit_finding'] = {'mad_threshold': args.threshold}
        if args.method:
            if 'transit_finding' not in config:
                config['transit_finding'] = {}
            config['transit_finding']['detrending_method'] = args.method

        catalogue = load_catalogue(args.catalogue, TEBC=args.tebc if hasattr(args, 'tebc') else False)
        finder = TransitFinder(
            catalogue=catalogue,
            sector_times=args.sector_times if hasattr(args, 'sector_times') else None,
            config=config
        )

        results = finder.process_directory(args.data_dir, output_file=args.output)
        print(f"Transit finding complete! Found {len(results)} events")

    elif args.command == 'compare-models':
        from mono_cbp import ModelComparator

        comparator = ModelComparator()
        results = comparator.compare_events(args.event_dir, output_file=args.output)
        print(f"Model comparison complete! Processed {len(results)} events")

    elif args.command == 'inject-retrieve':
        from mono_cbp import TransitInjector

        injector = TransitInjector(
            transit_models_path=args.models,
            catalogue=args.catalogue,
        )

        results = injector.run_injection_retrieval(
            args.data_dir,
            n_injections=args.n_injections,
            output_file=args.output
        )
        print(f"Injection-retrieval complete! Recovery rate: {injector.stats['recovery_rate']:.2%}")


if __name__ == '__main__':
    main()
