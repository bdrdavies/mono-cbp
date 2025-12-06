"""Command-line interface for mono-cbp."""

import argparse
import sys

def main():
    """Main CLI entry point.

    Provides command-line interface for running the mono-cbp pipeline and its components.
    """
    parser = argparse.ArgumentParser(
        description='mono-cbp: Pipeline for Detecting Transits of Circumbinary Planets in TESS Light Curves',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Run complete pipeline
    run_parser = subparsers.add_parser('run', help='Run complete pipeline')
    run_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV with eclipse and orbital parameters')
    run_parser.add_argument('--data-dir', required=True, help='Data directory')
    run_parser.add_argument('--sector-times', help='Path to sector times CSV')
    run_parser.add_argument('--output-dir', default='results', help='Output directory')
    run_parser.add_argument('--plot-dir', help='Directory to save vetting and Skye plots (if enabled in config)')
    run_parser.add_argument('--config', help='Path to configuration JSON file')
    run_parser.add_argument('--tebc', action='store_true', help='Use TEBC catalogue format (with *_2g and *_pf columns)')

    # Eclipse masking
    mask_parser = subparsers.add_parser('mask-eclipses', help='Mask eclipses only')
    mask_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV')
    mask_parser.add_argument('--data-dir', required=True, help='Data directory (files modified in-place)')
    mask_parser.add_argument('--config', help='Path to configuration JSON file')
    mask_parser.add_argument('--tebc', action='store_true', help='Use TEBC catalogue format (with *_2g and *_pf columns)')

    # Transit finding
    find_parser = subparsers.add_parser('find-transits', help='Find transits only')
    find_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV with eclipse and orbital parameters')
    find_parser.add_argument('--data-dir', required=True, help='Data directory')
    find_parser.add_argument('--sector-times', help='Path to sector times CSV')
    find_parser.add_argument('--output', default='transit_events.txt', help='Output file')
    find_parser.add_argument('--plot-dir', help='Directory to save vetting and Skye plots (if enabled in config)')
    find_parser.add_argument('--threshold', type=float, help='MAD threshold')
    find_parser.add_argument('--method', choices=['cb', 'cp'], help='Detrending method')
    find_parser.add_argument('--config', help='Path to configuration JSON file')
    find_parser.add_argument('--tebc', action='store_true', help='Use TEBC catalogue format (with *_2g and *_pf columns)')

    # Model comparison
    compare_parser = subparsers.add_parser('compare-models', help='Compare models for vetting')
    compare_parser.add_argument('--event-dir', required=True, help='Event snippets directory')
    compare_parser.add_argument('--output', default='classifications.csv', help='Output file')
    compare_parser.add_argument('--output-dir', help='Output directory (defaults to event-dir if not specified)')
    compare_parser.add_argument('--config', help='Path to configuration JSON file')

    # Injection-retrieval
    inject_parser = subparsers.add_parser('inject-retrieve', help='Run injection-retrieval')
    inject_parser.add_argument('--models', required=True, help='Path to transit models .npz')
    inject_parser.add_argument('--data-dir', required=True, help='Data directory')
    inject_parser.add_argument('--catalogue', required=True, help='Path to catalogue CSV with eclipse and orbital parameters')
    inject_parser.add_argument('--output', default='inj-ret_results.csv', help='Output file')
    inject_parser.add_argument('--output-dir', help='Output directory (defaults to data-dir if not specified)')
    inject_parser.add_argument('--n-injections', type=int, default=None, help='Number of injections per model (defaults to config value)')
    inject_parser.add_argument('--config', help='Path to configuration JSON file')
    inject_parser.add_argument('--tebc', action='store_true', help='Use TEBC catalogue format (with *_2g and *_pf columns)')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Import here to avoid slow startup
    if args.command == 'run':
        from mono_cbp import MonoCBPPipeline
        import json

        print("\n=== Running Complete Pipeline ===")
        print(f"Catalogue: {args.catalogue}")
        print(f"Data directory: {args.data_dir}")
        print(f"Output directory: {args.output_dir if hasattr(args, 'output_dir') else './results'}")
        if hasattr(args, 'plot_dir') and args.plot_dir:
            print(f"Plot directory: {args.plot_dir}")

        config = {}
        if hasattr(args, 'config') and args.config:
            print(f"Configuration: {args.config}")
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

        print("\nStarting pipeline...\n")

        # Prepare kwargs for pipeline steps
        pipeline_kwargs = {}
        if hasattr(args, 'plot_dir') and args.plot_dir:
            pipeline_kwargs['find_transits_kwargs'] = {'plot_output_dir': args.plot_dir}

        results = pipeline.run(**pipeline_kwargs)

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
        import json

        print("\n=== Masking Eclipses ===")
        print(f"Catalogue: {args.catalogue}")
        print(f"Data directory: {args.data_dir}")
        print("Note: Files will be modified in-place\n")

        # Load config from file if provided (for potential future use)
        config = {}
        if hasattr(args, 'config') and args.config:
            print(f"Configuration: {args.config}")
            with open(args.config) as f:
                config = json.load(f)

        catalogue = load_catalogue(args.catalogue, TEBC=args.tebc if hasattr(args, 'tebc') else False)
        masker = EclipseMasker(catalogue, data_dir=args.data_dir)

        print("Starting eclipse masking...\n")
        masker.mask_all()
        print("\nEclipse masking complete!")

    elif args.command == 'find-transits':
        from mono_cbp import TransitFinder
        from mono_cbp.utils import load_catalogue
        import json

        print("\n=== Finding Transits ===")
        print(f"Catalogue: {args.catalogue}")
        print(f"Data directory: {args.data_dir}")
        print(f"Output file: {args.output}")
        if hasattr(args, 'plot_dir') and args.plot_dir:
            print(f"Plot directory: {args.plot_dir}")
        if args.threshold:
            print(f"MAD threshold: {args.threshold}")
        if args.method:
            print(f"Detrending method: {args.method}")

        # Load config from file if provided
        config = {}
        if hasattr(args, 'config') and args.config:
            print(f"Configuration: {args.config}")
            with open(args.config) as f:
                config = json.load(f)

        # Override config with command-line arguments
        if args.threshold:
            if 'transit_finding' not in config:
                config['transit_finding'] = {}
            config['transit_finding']['mad_threshold'] = args.threshold
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

        print("\nStarting transit finding...\n")
        results = finder.process_directory(
            args.data_dir,
            output_file=args.output,
            plot_output_dir=args.plot_dir if hasattr(args, 'plot_dir') and args.plot_dir else None
        )
        print(f"\nTransit finding complete! Found {len(results)} events")

    elif args.command == 'compare-models':
        from mono_cbp import ModelComparator
        import json

        print("\n=== Comparing Models for Vetting ===")
        print(f"Event directory: {args.event_dir}")
        print(f"Output file: {args.output}")
        if hasattr(args, 'output_dir') and args.output_dir:
            print(f"Output directory: {args.output_dir}")
        else:
            print(f"Output directory: {args.event_dir} (default)")

        # Load config from file if provided
        config = {}
        if hasattr(args, 'config') and args.config:
            print(f"Configuration: {args.config}")
            with open(args.config) as f:
                config = json.load(f)

        comparator = ModelComparator(config=config if config else None)

        print("\nStarting model comparison...\n")
        results = comparator.compare_events(
            args.event_dir,
            output_file=args.output,
            output_dir=args.output_dir if hasattr(args, 'output_dir') and args.output_dir else None
        )
        print(f"\nModel comparison complete! Processed {len(results)} events")

    elif args.command == 'inject-retrieve':
        from mono_cbp import TransitInjector
        import json

        # Load config from file if provided
        config = {}
        if hasattr(args, 'config') and args.config:
            print(f"Configuration: {args.config}")
            with open(args.config) as f:
                config = json.load(f)

        # Determine n_injections (from args or config)
        from mono_cbp.config import get_default_config, merge_config
        merged_config = merge_config(config, get_default_config()) if config else get_default_config()
        n_injections_display = args.n_injections if args.n_injections is not None else merged_config['injection_retrieval']['n_injections']

        print("\n=== Running Injection-Retrieval ===")
        print(f"Transit models: {args.models}")
        print(f"Catalogue: {args.catalogue}")
        print(f"Data directory: {args.data_dir}")
        print(f"Number of injections per model: {n_injections_display}{' (from config)' if args.n_injections is None else ''}")
        print(f"Output file: {args.output}")
        if hasattr(args, 'output_dir') and args.output_dir:
            print(f"Output directory: {args.output_dir}")
        else:
            print(f"Output directory: {args.data_dir} (default)")


        injector = TransitInjector(
            transit_models_path=args.models,
            catalogue=args.catalogue,
            config=config if config else None,
            TEBC=args.tebc if hasattr(args, 'tebc') else False
        )

        print("\nStarting injection-retrieval...\n")
        results = injector.run_injection_retrieval(
            args.data_dir,
            n_injections=args.n_injections,
            output_file=args.output,
            output_dir=args.output_dir if hasattr(args, 'output_dir') and args.output_dir else None
        )
        print(f"\nInjection-retrieval complete! Processed {len(results)} injections")


if __name__ == '__main__':
    main()